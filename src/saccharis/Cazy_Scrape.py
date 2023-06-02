###############################################################################
#
#
# SACCHARIS 2.0 author: Alexander Fraser, https://github.com/AlexSCFraser
# Original author for SACCHARIS 1.0: Dallas Thomas
# License: GPL v3
###############################################################################
import math
import time
from enum import Enum
import re
import os
import csv
import json
# ignore warnings that these parsers aren't used by IDE, they ARE used
# noinspection PyUnresolvedReferences
from html.parser import HTMLParser
# noinspection PyUnresolvedReferences
import lxml
import requests
from bs4 import BeautifulSoup

from saccharis.utils.PipelineErrors import PipelineException
from saccharis.utils.AdvancedConfig import load_from_env

NCBI_DEFAULT_DELAY = 0.3  # this is a delay time for ncbi queries. Without it, results may be incomplete
# count ncbi exceptions, so we can terminate if too many failures occur. Too many failures probably means NCBI is down.
ncbi_exception_count = 0
NCBI_EXCEPTION_MAX_TRIES = 100


class Mode(Enum):
    CHARACTERIZED = 0
    ALL_CAZYMES = 1
    STRUCTURE = 2


class Domain(Enum):
    # Using these values for a bitmask comparison to easily combine groups, which is why we use hex values
    # ALL = 0b11111 # correct value for all, but is omitted so that list comprehension works correctly
    ARCHAEA = 0b00001
    BACTERIA = 0b00010
    EUKARYOTA = 0b00100
    VIRUSES = 0b01000
    UNCLASSIFIED = 0b10000


class NCBIException(PipelineException):
    def __init__(self, msg):
        super().__init__(msg)


class CazyException(PipelineException):
    def __init__(self, msg):
        super().__init__(msg)


class HTMLGetter:

    def __init__(self):
        self.fix_spaces = re.compile("> <")
        self.fix_white = re.compile("#FFFFFF")

    def get_clean_html_text(self, url_cazy, tries=0):
        try:
            r = requests.get(url_cazy)
        except requests.ConnectionError as error:
            raise PipelineException("Connection error, cazy.org might be down right now.") from error
        except requests.RequestException as error:
            raise PipelineException("HTTP request error, cazy.org might be down or overloaded right now.") from error

        if r.status_code != 200 and tries < 3:
            print(f"Bad http response {r.status_code} from CAZy.org, retrying...")
            time.sleep(3)
            return self.get_clean_html_text(url_cazy, tries+1)
        elif r.status_code != 200:
            print(f"ERROR: HTTP response status code {r.status_code}")
            print("Max tries to retrieve data from CAZy.org reached.")
            print(f"CAZy.org might be down, or failing to fulfill requests to {url_cazy}")
            print("Aborting current pipeline iteration...")
            raise PipelineException(f"Bad HTTP response code {r.status_code} from {url_cazy}, max tries exceeded.")

        # remove bad tag spacing and other tag fixes
        clean_text = self.fix_spaces.sub("><", r.text)
        clean_text = self.fix_white.sub("#ffffff", clean_text)

        return clean_text


def valid_genbank(string_to_check, verbose=False):
    if string_to_check is None or string_to_check == "":
        return False

    # This matches several standard genbank accession formats, e.g.
    #  Standard Accession formats:
    #  [one-letter alphabetical prefix][five digits][.][version number]
    #  [two-letter alphabetical prefix][six digits][.][version number]
    #  [three-letter alphabetical prefix][5 digits][.][version number]
    #  Reference:   https://support.nlm.nih.gov/knowledgebase/article/KA-03436/en-us

    #  Official documentation on what constitutes a valid accession seems out of date. I expect there may be longer
    #  alphabet character prefixes in the future(e.g. 4+ letters), so I just catch all alphanumeric characters with
    #  the regex.

    #  I also catch an additional accession format for non-redundant sequences, again loosely matching the correct
    #  number of alphanumeric characters only since the format may change in the future:
    #  [two letter alphabetical prefix][_][9 digits][.][version number]
    #  Example: WP_010248927.1
    #  Reference:   https://www.ncbi.nlm.nih.gov/refseq/about/nonredundantproteins/

    #  Note that this function does NOT return true for the expanded format, used for Whole Genome Shotgun (WGS),
    #  Transcriptome Shotgun Assembly (TSA), and Targeted Locus Study (TLS) sequencing projects. Those use a six-letter
    #  Project Code prefix and a two-digit Assembly-Version number followed by 7, 8, or 9 digits
    #  (for example, AAAAAA020000001).
    #  https://www.nlm.nih.gov/pubs/techbull/so18/brief/so18_genbank_expanded_accession_formats.html

    if re.match(r"\w{8}|\w{6}|\w{12}\.\d+", string_to_check):
        return True

    if verbose:
        print(f"\"{string_to_check}\" is not a valid genbank accession identifier. If this is in error(i.e. this string"
              f" is in fact a valid genbank accession), please report this as a bug to the developer via the github "
              f"page and/or email.")
    return False


def format_list(accession_list):
    genbank_list = ""
    for accession_id in accession_list:
        if accession_id is not None:
            genbank_list += (accession_id + ',')
    if len(genbank_list) > 0:
        genbank_list = genbank_list[0:-1]   # remove ',' at end
    return genbank_list


def ncbi_query(accession_list, api_key, ncbi_email, ncbi_tool, verbose=False):
    genbank_list = format_list(accession_list)

    # Set up the Query URL
    # Set up a search out to the eSearch program:
    #    - db is protein and search term is left blank for now <term>
    email_string = '&email=' + ncbi_email if ncbi_email else ""
    tool_string = '&tool=' + ncbi_tool if ncbi_tool else ""
    api_key_string = f"&api_key={api_key}" if api_key else ""
    # todo: consider checking for valid API_key here
    utils = "http://www.ncbi.nlm.nih.gov/entrez/eutils"
    # TODO: add email and tool fields to ncbi requests
    base_url = utils + "/esearch.fcgi?db=protein" + api_key_string
    esearch = base_url + '&term='

    # Submit the search to retrieve a count of total number of sequences
    try:
        time.sleep(NCBI_DEFAULT_DELAY)
        esearch_result1 = requests.get(esearch + genbank_list)
    # todo: consider catching specific exceptions here. These are intermittent and not repeatable, since they happen
    #  when the NCBI server has errors, so I am not sure which specific exceptions to catch.
    except Exception as e:
        print(e.args[0])
        msg = "Error querying NCBI. NCBI might be down, try again later.\nFailed NCBI request #1.\n"
        raise NCBIException(msg) from e

    # Extract the count and submit search again to retrieve XML based results
    # - set the number of results we want to count <retmax>
    valid_accession_count = len(accession_list)
    result1 = BeautifulSoup(esearch_result1.text, features='xml')

    # Remove accession numbers that were not found, count valid, rebuild the list for querying
    bad_accessions = result1.find_all('PhraseNotFound')
    for item in bad_accessions:
        if verbose:
            print("\nWARNING: NCBI DATA MISSING")
            print("Genbank accession:", item.text, "\n")
        accession_list.remove(item.text)
        valid_accession_count -= 1

    # # Note: The counting sum below does not always work. Sometimes there are substantially fewer <Count></Count>
    # #       tags than there are valid entries. I have no idea why this is, but for now we just set max_ret to the
    # #       valid_accession_count. This may fail if there are multiple entries off a single protein.
    # #       NOTE: The number in between count tags should be the total valid entries (for queries without errors)
    # #             or a series of count tags that are either 0(entry not found) or 1(entry found). There should
    # #             be no situations where multiple FASTA sequences are returned from a single entry in a query for
    # #             accessions numbers kept in the CAZy database.
    # #
    # count_result1 = result1.find_all('Count')
    # max_ret = 0
    # for i, tag in enumerate(count_result1):
    #     # if int(tag.contents[0]) == 0:
    #     #     accession_list[i] = None
    #     #     valid_accession_count -= 1
    #     max_ret += int(tag.contents[0])    # not sure if sum is what we want or not?

    max_ret = valid_accession_count
    genbank_list = format_list(accession_list)

    esearch = base_url + '&retmax=' + str(max_ret) + '&term='
    try:
        time.sleep(NCBI_DEFAULT_DELAY)
        esearch_result2 = requests.get(esearch + genbank_list + '&usehistory=y')
    # todo: consider catching specific exceptions here. These are intermittent and not repeatable, since they happen
    #  when the NCBI server has errors, so I have no idea which specific exceptions to catch.
    except Exception as e:
        print(e.args[0])
        raise NCBIException("Error querying NCBI. NCBI might be down, try again later.\nFailed NCBI request #2")

    result2 = BeautifulSoup(esearch_result2.text, features='xml')
    if result2.find('QueryKey') is None and result2.find('querykey') is None:
        raise NCBIException("ERROR: NCBI query Key not found. Usually this means query size is too large.")
    if result2.find('QueryKey') is None:
        query_key = result2.find('querykey').text
    else:
        query_key = result2.find('QueryKey').text
    if result2.find('WebEnv') is None:
        web_env = result2.find('webenv').text
    else:
        web_env = result2.find('WebEnv').text

    # Fetch the Fasta data from NCBI using the esearch results
    # $efetch = $utils . '/efetch.fcgi?db=protein&query_key=' . $key . '&WebEnv='
    #             . $web . '&rettype=fasta&retmode=text';
    # base_url = utils + '/efetch.fcgi?db=protein&email=' + ncbi_email + '&tool=' + ncbi_tool + '&api_key='+api_key
    # TODO: add email and tool fields to ncbi requests
    base_url = utils + '/efetch.fcgi?db=protein' + api_key_string
    efetch = base_url + '&query_key=' + query_key + '&WebEnv=' + web_env + '&rettype=fasta&retmode=text'

    try:
        time.sleep(NCBI_DEFAULT_DELAY)
        efetch_result = requests.get(efetch)
        result_count = efetch_result.text.count(">")
    # todo: consider catching specific exceptions here. These are intermittent and not repeatable, since they happen
    #  when the NCBI server has errors, so I have no idea which specific exceptions to catch.
    except Exception as e:
        print(e.args[0])
        raise NCBIException("HTTP error querying NCBI. NCBI might be down, try again later.\nFailed NCBI request #3")

    if valid_accession_count != result_count:
        raise NCBIException(f"Incomplete NCBI query FASTA results: {result_count}/{valid_accession_count} returned")

    # Returns empty result if fetch failed
    if efetch_result.text.__contains__("<ERROR>Empty result - nothing to do</ERROR>"):
        print("ERROR: NCBI Fetch failed.")
        return "", 0

    # Remove double newline between each of the sequences
    fasta_out = re.sub("\n+", "\n", efetch_result.text)

    # Replaces weird NCBI accessions like 'sp|<ACCESSION>|', 'prf|<ACCESSION>|', or 'pir||<ACCESSION>' with <ACCESSION>
    # This retains the '>' at the beginning of text lines in the FASTA data and retains auxillary data after the '|' by
    # concatenating it with space separations.
    # Another way to think about this is there is a pattern of ">(IDENTIFIER)|(ACCESSION)|(AUXILLARY_DATA)" and we
    # delete the identifier and first '|' character while replacing subsequent '|' characters in the first word with
    # spaces.
    if fasta_out.__contains__('|'):
        lines = fasta_out.split('\n')
        for i, row in enumerate(lines):
            words = row.split(' ')
            if words[0].__contains__('|'):
                words[0] = re.sub("\|+", "|", words[0])
                accession_array = words[0].split('|')
                accession_array.pop(0)
                words[0] = ' '.join(accession_array)
                lines[i] = '>' + ' '.join(words)
        fasta_fixed = '\n'.join(lines)
        fasta_out = fasta_fixed

    if fasta_out.__contains__('|'):
        print(f"WARNING: Probable parsing error on accession containing a | character, "
              f"Please report this as a bug to the developer/maintainer through github.")

    return fasta_out, result_count


def cazy_query(family, cazy_folder, scrape_mode, get_fragments, verbose, domain_mode):
    print("Downloading", family, "Data from CAZy database...")

    url_cazy = "http://www.cazy.org/"+family+"_" + \
               ("structure" if scrape_mode == Mode.STRUCTURE or family.__contains__("CBM") else "characterized") \
               + ".html"
    html_get = HTMLGetter()
    clean_text = html_get.get_clean_html_text(url_cazy)

    # get number of CAZymes in family
    soup = BeautifulSoup(clean_text, "html.parser")
    try:
        count_string = soup.find("span", id="line_actif").contents[1].strip()  # span with this id contains the count
        if count_string.__contains__("-"):  # - is in this string when getting from structure tab, "cryst" in tabname
            count_string = count_string.split("-")[0].strip()
            count_string = count_string[1:].strip()
        else:
            count_string = count_string[1:-1]
    except AttributeError as soup_error:
        print("EXCEPTION/DEBUG MESSAGE:", soup_error.args[0])
        raise CazyException(f"ERROR: No entries for family {family} in the CAZy html response. CAZy query of "
                            f"{family} aborted. If there ARE entries for this family, please report this as a bug.") \
            from soup_error

    count = int(count_string)
    counted = 0
    cazy_retrieved = 0
    cazy_duplicate = 0
    cazy_fragments = 0
    cazy_missing = 0
    cazy_added = 0

    cazymes = {}
    genbank_duplicates = []

    # loop through all pages of characterized CAZymes for selected family
    while counted < count:
        if counted > 0:
            # load new page on second iteration and above
            url = url_cazy + f"?debut_FUNC={str(counted)}#pagination_FUNC"
            clean_text = html_get.get_clean_html_text(url)
            soup = BeautifulSoup(clean_text, "html.parser")

        # find and filter table entries
        # tables = soup.find_all('table')
        # # TODO: improve this filter? search <tr> children for e.g. genbank accession pattern or EC #?
        # rows = tables[1].find_all("tr", attrs={'bgcolor': "#ffffff"})       # filter only rows that are CAZyme entries
        rows = soup.find_all("tr",  attrs={'bgcolor': "#ffffff"})
        header_row = soup.select("tr[id='line_titre'] > td:not([colspan])")
        col_idx = {"Protein Name": None,
                   "EC#": None,
                   "Reference": None,
                   "Organism": None,
                   "GenBank": None,
                   "Uniprot": None,
                   "PDB/3D": None,
                   "Subf": None}
        for i, col in enumerate(header_row):
            col_name = col.text.strip()
            if col_name.__contains__(" Carbohydrate Ligands"):
                col_name = col_name.split(' ')[0]

            if col_name and not col_name.__contains__("Å") and col_idx[col_name] is None:
                col_idx[col_name] = 2*i + 1

        for entry in rows:
            cazy_retrieved += 1
            # GenBank
            genbank = None
            for child in entry.contents[col_idx["GenBank"]]:
                if genbank is None and child.name != "br":
                    genbank = child.text.strip()
                elif child.name != "br":
                    genbank_duplicates.append(child.text.strip())
            # todo: refactor the rest of the attributes to syntax like above for clarity
            if scrape_mode == Mode.STRUCTURE or family.__contains__("CBM"):  # CBM families have no characterized page
                protein_name = entry.contents[1].contents[0].strip()                    # protein name
                ec_num = ""
                try:
                    for ec_blob in entry.contents[3]:                                   # iterate through all ec nums
                        if len(ec_blob.contents) > 0:
                            # ec_num += ec_blob.contents[0].attrs['href'].strip()       # concat ec_num url
                            ec_num += ec_blob.contents[0].contents[0].strip()           # concat ec_num
                        else:
                            ec_num += " "                                               # concat delimiter for ec_num
                except AttributeError:
                    ec_num = None
                ref_art = None        # reference article - ITERATE THIS?

                try:
                    org_name = entry.contents[5].contents[1].contents[0].contents[0].strip()    # Organism
                except AttributeError:
                    org_name = entry.contents[5].contents[1].contents[0].strip()                # Organism
                try:
                    # uniprot = entry.contents[11].contents[0].text.strip()
                    uniprot = ""
                    for uni_blob in entry.contents[9]:
                        text = uni_blob.text.strip()
                        if len(text) > 0:
                            uniprot += text
                        else:
                            uniprot += " "
                    uniprot = uniprot.strip()
                    if uniprot == "":
                        uniprot = None
                except AttributeError:
                    uniprot = None
                try:
                    pdb = ""
                    for pdb_blob in entry.contents[11].contents[1]:
                        if pdb_blob.text == "\n":
                            pdb += " "
                        else:
                            text = pdb_blob.contents[1].text.strip()
                            if len(text) > 0:
                                pdb += text
                            else:
                                pdb += " "
                    pdb = pdb.strip()
                    if pdb == "":
                        pdb = None
                except AttributeError:
                    pdb = None          # todo: check if there is another nested tag to check and get id not url
                try:
                    subfamily = entry.contents[13].contents[0].strip()  # subf
                except AttributeError:
                    subfamily = None
                except IndexError:
                    subfamily = None
            else:   # scrape mode == ALL or CHARACTERIZED
                protein_name = entry.contents[1].contents[0].strip()                    # protein name
                ec_num = ""
                for ec_blob in entry.contents[3]:                                       # iterate through all ec nums
                    if len(ec_blob.contents) > 0:
                        # ec_num += ec_blob.contents[0].attrs['href'].strip()                   # concat ec_num url
                        ec_num += ec_blob.contents[0].contents[0].strip()                       # concat ec_num
                    else:
                        ec_num += " "                                                   # concat delimiter for ec_num
                try:
                    ref_art = entry.contents[5].contents[0].attrs['href'].strip()    # reference article - ITERATE THIS?
                except AttributeError:
                    ref_art = None
                try:
                    org_name = entry.contents[7].contents[1].contents[0].contents[0].strip()    # Organism
                except AttributeError:
                    org_name = entry.contents[7].contents[1].contents[0].strip()                # Organism
            # print(entry.contents[11].contents[0].attrs['href'])                   # Uniprot # Add code to handle blank
            # print(entry.contents[13].contents[0])                                 # PDB #add code to handle blank
            # print(tables[1].contents[7].contents[15].contents[0])                 # subf
                try:
                    subfamily = entry.contents[15].contents[0].strip()  # subf
                except AttributeError:  # Attribute error does not normally occur, might remove for clarity?
                    subfamily = None
                except IndexError:
                    subfamily = None
                try:
                    # uniprot = entry.contents[11].contents[0].text.strip()
                    uniprot = ""
                    for uni_blob in entry.contents[11]:
                        text = uni_blob.text.strip()
                        if len(text) > 0:
                            uniprot += text
                        else:
                            uniprot += " "
                    uniprot = uniprot.strip()

                    if uniprot == "":
                        uniprot = None
                except AttributeError:
                    uniprot = None
                try:
                    pdb = ""
                    for pdb_blob in entry.contents[13]:
                        text = pdb_blob.text.strip()
                        if len(text) > 0:
                            pdb += text
                        else:
                            pdb += " "
                    pdb = pdb.strip()
                    if pdb == "":
                        pdb = None
                except AttributeError:
                    pdb = None

            # print(protein_name)
            # print(ec_num)
            # print(ref_art)
            # print(org_name)
            # print(genbank)
            # print("\n")

            if valid_genbank(genbank, verbose) and genbank not in cazymes and \
                    (get_fragments or not protein_name.__contains__("fragment")):
                # todo: change cazymes from a dict of lists to a dict of dicts (or a dict of Namespace objects? or dict
                #  of custom class?) to get named json categories in output. THIS WILL BREAK DATA IMPORT INTO R SCRIPT
                cazymes[genbank] = [protein_name, ec_num, org_name, None, uniprot, pdb, subfamily]  # None is domain, filled later
                cazy_added += 1
            else:
                if genbank is not None and genbank != '' and genbank in cazymes:
                    cazy_duplicate += 1
                    skip_msg = "DUPLICATE - CHARACTERIZED CAZYME NOT ADDED:"
                elif genbank is not None and protein_name.__contains__("fragment"):
                    cazy_fragments += 1
                    skip_msg = "FRAGMENT - CHARACTERIZED CAZYME NOT ADDED:"
                else:   # else we assume genbank data was missing
                    # todo: If genbank is None and we are in structure mode, should we try to get sequence from pdb?
                    cazy_missing += 1
                    skip_msg = "MISSING - CHARACTERIZED CAZYME NOT ADDED:"
                if verbose:
                    print("\n" + skip_msg)
                    print("Protein name:", protein_name)
                    print("EC #:", ec_num)
                    print("Article:", ref_art)
                    print("Organism:", org_name)
                    print("Genbank:", genbank, "\n")

        counted += 100

    total_count = 0
    uncharacterized_added = 0
    uncharacterized_duplicate = 0
    unchar_char_duplicate = 0

    # download "all" file, used getting uncharacterized and also identifying domain
    url_all = "http://www.cazy.org/IMG/cazy_data/" + family + ".txt"
    full_list = requests.get(url_all)
    if full_list.status_code != 200:
        raise PipelineException(f"Bad HTTP response code {full_list.status_code} from {url_all}")
    full_path = os.path.join(cazy_folder, family + '_full_list.txt')
    with open(full_path, 'wb') as listfile:
        listfile.write(full_list.content)

    if scrape_mode == Mode.ALL_CAZYMES:
        all_cazymes = cazymes
        with open(full_path, 'r', newline='\n') as csvfile:
            entry_reader = csv.reader(csvfile, delimiter='\t')
            for row in entry_reader:
                total_count += 1
                cazyme_class = row[0]
                domain = row[1]
                organism = row[2]
                genbank = row[3]
                if genbank not in all_cazymes and genbank not in genbank_duplicates and genbank is not None and genbank != "":
                    uncharacterized_added += 1

                    # genbank_query = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/" + )
                    all_cazymes[genbank] = [f"Uncharacterized {cazyme_class}", None, organism, domain, None, None]
                # we check genbank not in cazymes to prevent reporting characterized as duplicates
                elif genbank not in cazymes:
                    uncharacterized_duplicate += 1
                    if verbose:
                        print("DUPLICATE - UNCHARACTERIZED CAZYME NOT ADDED:")
                        print("Organism:", organism)
                        print("Genbank:", genbank)
                        print("File Line #:", entry_reader.line_num)
                        print("\n")
                elif genbank in cazymes:
                    # add domain to characterized entry
                    unchar_char_duplicate += 1
                    cazymes[genbank][3] = domain
                else:
                    print("ERROR: Uncharacterized entry not parsed correctly, please report this as a bug.")
        cazymes = all_cazymes
        uncharacterized_retrieved = total_count-unchar_char_duplicate
    else:
        uncharacterized_retrieved = 0
        # add domain to characterized entries
        with open(full_path, 'r', newline='\n') as csvfile:
            entry_reader = csv.reader(csvfile, delimiter='\t')
            for row in entry_reader:
                total_count += 1
                # cazyme_class = row[0]
                domain = row[1]
                # organism = row[2]
                genbank = row[3]
                if genbank in cazymes:
                    cazymes[genbank][3] = domain

    # Filter for correct domain
    wrong_domain_characterized = 0
    wrong_domain_uncharacterized = 0
    domain_cazymes = {}
    for item in cazymes:
        # bitwise comparison against bitmask
        if cazymes[item][3] and Domain[cazymes[item][3].upper()].value & domain_mode:
            domain_cazymes[item] = cazymes[item]
        else:  # invalid domain
            if cazymes[item][1] or cazymes[item][4] or cazymes[item][5]:
                wrong_domain_characterized += 1
            else:
                wrong_domain_uncharacterized += 1

    # # This list comprehension works fine, but I changed to the for loop above to track wrong_domain vars
    # domain_cazymes = {item for item in cazymes if Domain[cazymes[item][3].upper()].value & domain_mode}
    # wrong_domain = len(cazymes) - len(domain_cazymes)

    cazy_added -= wrong_domain_characterized
    uncharacterized_added -= wrong_domain_uncharacterized

    print("Done retrieving cazymes data from CAZy!")
    # stats = [cazy_retrieved, cazy_added, cazy_duplicate, cazy_fragments, cazy_missing,
    #          uncharacterized_retrieved, uncharacterized_added, uncharacterized_duplicate]
    stats = [cazy_retrieved, cazy_added, cazy_duplicate, cazy_fragments, cazy_missing, wrong_domain_characterized,
             uncharacterized_retrieved, uncharacterized_added, uncharacterized_duplicate, wrong_domain_uncharacterized]

    if stats[0] != stats[1] + stats[2] + stats[3] + stats[4] + stats[5] or stats[6] != stats[7] + stats[8] + stats[9]:
        print("WARNING / ERROR - Summary statistics on CAZy retrieval do not add up, the statistics are not reliable, "
              "please file a bug report with the developer.")

    return domain_cazymes, stats


def main(family, cazy_folder, scrape_mode, get_fragments=False, verbose=False, force_update=False, ncbi_query_size=200,
         domain_mode=0b11111, gui_active=False):
    api_key, ncbi_email, ncbi_tool = load_from_env(skip_ask=gui_active)
    # Folder and file output setup
    fasta_file = os.path.join(cazy_folder, f"{family}_{scrape_mode.name}_cazy.fasta")
    data_file = os.path.join(cazy_folder, f"{family}_{scrape_mode.name}_data.json")
    stats_file = os.path.join(cazy_folder, f"{family}_{scrape_mode.name}_stats.json")

    if os.path.isfile(fasta_file) and not force_update:
        print("Loading data from previous run. \nIf you wish to run a new query to the CAZy database, run SACCHARIS 2"
              " with --fresh")
        with open(data_file, 'r', encoding='utf-8') as f:
            cazymes = json.loads(f.read())
        with open(stats_file, 'r', encoding='utf-8') as f:
            stats = json.loads(f.read())
        return fasta_file, cazymes, stats

    cazymes, cazy_stats = cazy_query(family, cazy_folder, scrape_mode, get_fragments, verbose, domain_mode)

    # Take the accession numbers from the dict, convert to list, and query genbank in batches
    accession_list = list(cazymes.keys())
    accession_count = len(accession_list)
    queried = 0
    retrieved = 0
    fasta_data = ""
    print("Querying NCBI...", end='')
    while queried < accession_count:
        # fasta_output, success_count = ncbi_query(accession_list[queried:queried+ncbi_query_size])
        try:
            fasta_output, success_count = ncbi_query(accession_list[queried:queried+ncbi_query_size], api_key,
                                                     ncbi_email, ncbi_tool)
            fasta_data += fasta_output
            retrieved += success_count
            queried += ncbi_query_size
            queried = min(queried, accession_count)
            print("\rQuerying NCBI: %d/%d entries processed..." % (queried, len(accession_list)), end='')
        except NCBIException as error:
            global ncbi_exception_count
            ncbi_exception_count += 1
            print("WARNING: MISSING FASTA DATA FROM NCBI")
            print(error.msg)
            if ncbi_exception_count < NCBI_EXCEPTION_MAX_TRIES:
                ncbi_query_size = math.ceil(ncbi_query_size/2)
                print(f"INFO: Automatically reducing query size to {ncbi_query_size} and retrying...\n")
            else:
                raise PipelineException("Exceeded maximum failed NCBI query attempts. This probably means NCBI servers "
                                        "are down or not responding. Consider waiting a while and trying again.\n"
                                        "\tIf this problem persists, please contact the developer to investigate.\n")
    print("\rQuerying NCBI...Done!                                  ")  # whitespace to overwrite "entries processed"
    cazy_stats.append(queried)
    cazy_stats.append(retrieved)

    # write standard fasta file
    with open(fasta_file, 'w') as f:
        f.write(fasta_data)
    # write data file of all the ancillary data as a dict
    with open(data_file, 'w', encoding='utf-8') as f:
        json.dump(cazymes, f, ensure_ascii=False, indent=4)
    # write stats file
    with open(stats_file, 'w', encoding='utf-8') as f:
        json.dump(cazy_stats, f, ensure_ascii=False, indent=4)

    print("Characterized entries retrieved from CAZy database:", cazy_stats[0])
    print("Characterized entries added to dataset:", cazy_stats[1])
    if cazy_stats[2] > 0:
        print("Characterized entries with duplicate accessions not added:", cazy_stats[2])
    if cazy_stats[3] > 0:
        print("Characterized entries that are fragments not added:", cazy_stats[3])
    if cazy_stats[4] > 0:
        print("Characterized entries with missing accession not added:", cazy_stats[4])
    if cazy_stats[5] > 0:
        print("Characterized entries with wrong domain not added:", cazy_stats[5])
    if scrape_mode == Mode.ALL_CAZYMES:
        print("Uncharacterized entries retrieved from CAZy database:", cazy_stats[6])
        print("Uncharacterized entries added to dataset:", cazy_stats[7])
        if cazy_stats[7] > 0:
            print("Uncharacterized entries with duplicate accessions not added:", cazy_stats[8])
        if cazy_stats[8] > 0:
            print("Uncharacterized entries with wrong domain not added:", cazy_stats[9])
    print("Total number of queried accessions:", cazy_stats[10])
    print("Total number of results from NCBI:", cazy_stats[11])
    if cazy_stats[10] != cazy_stats[11]:
        print("WARNING: MISSING DATA FROM NCBI. SOME FASTA RESULTS COULD NOT BE QUERIED!!!")
    if cazy_stats[2] > 0 or cazy_stats[3] > 0 or cazy_stats[4] > 0 or cazy_stats[7] > 0 \
            or cazy_stats[8] != cazy_stats[9]:
        print("Details on duplicate/fragment/missing accessions printed above." if verbose else
              " Enable verbose argument for more details on duplicate/fragment/missing accessions.")
    print("\n")

    return fasta_file, cazymes, cazy_stats


if __name__ == "__main__":
    # test args
    # test_family = "PL9"
    # test_family = "CBM4"
    test_family = "GH156"
    test_scrape_mode = Mode.CHARACTERIZED
    test_get_fragments = False
    test_verbose = True
    output_folder = os.path.join(os.getcwd(), "output")
    family_folder = os.path.join(output_folder, test_family)
    test_group_folder = os.path.join(output_folder, test_family, test_scrape_mode.name)

    if not os.path.isdir(family_folder):
        os.mkdir(family_folder, 0o755)
    if not os.path.isdir(test_group_folder):
        os.mkdir(test_group_folder, 0o755)

    main(test_family, test_group_folder, test_scrape_mode, test_get_fragments, test_verbose, force_update=True)
