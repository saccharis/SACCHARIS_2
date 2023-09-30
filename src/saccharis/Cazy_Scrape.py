###############################################################################
#
#
# SACCHARIS 2.0 author: Alexander Fraser, https://github.com/AlexSCFraser
# Original author for SACCHARIS 1.0: Dallas Thomas
# License: GPL v3
###############################################################################
import time
from enum import Enum
import re
import os
import csv
import json
from io import StringIO
from json import JSONDecodeError
from logging import Logger, getLogger
# ignore warnings that these parsers aren't used by IDE, they ARE used
# noinspection PyUnresolvedReferences
from html.parser import HTMLParser
# noinspection PyUnresolvedReferences
import lxml
import requests
from Bio import SeqIO
from bs4 import BeautifulSoup

from saccharis.NCBIQueries import valid_genbank_gene, ncbi_protein_query
from saccharis.utils.UserInput import ask_yes_no
from saccharis.utils.PipelineErrors import PipelineException
from saccharis.utils.AdvancedConfig import load_from_env
from saccharis.utils.PipelineErrors import CazyException
from saccharis.utils.Formatting import CazymeMetadataRecord


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


class HTMLGetter:

    def __init__(self):
        self.fix_spaces = re.compile("> <")
        self.fix_white = re.compile("#FFFFFF")

    def get_clean_html_text(self, url_cazy: str, tries: int = 0, logger: Logger = getLogger()):
        try:
            r = requests.get(url_cazy)
        except requests.ConnectionError as error:
            raise PipelineException("Connection error, cazy.org might be down right now.") from error
        except requests.RequestException as error:
            raise PipelineException("HTTP request error, cazy.org might be down or overloaded right now.") from error

        if r.status_code == 503 and tries < 5:
            # service unavailable, check for retry-after header and retry
            try:
                timeout = int(r.headers["retry-after"]) + 1
                msg = f"CAZy.org is unavailable, expected to be back up in {timeout} seconds..."
                logger.warning(msg)
            except KeyError as err:
                timeout = 60*tries
                msg = f"CAZy.org is unavailable, no expected time when the website is back up returned. \n" \
                      f"Retrying in {timeout} seconds..."

            logger.warning(msg)
            time.sleep(timeout)
            return self.get_clean_html_text(url_cazy, tries + 1)
        if r.status_code == 429 and tries < 5:
            # too many requests, check for retry-after header and retry
            try:
                timeout = int(r.headers["retry-after"]) + 1
                msg = f"Too many requests to CAZy.org, must wait {timeout} seconds before retrying..."
            except KeyError:
                timeout = 30*tries
                msg = f"Too many requests to CAZy.org, no timeout present, waiting {timeout} seconds..."

            logger.warning(msg)
            time.sleep(timeout)
            return self.get_clean_html_text(url_cazy, tries + 1)
        if r.status_code != 200 and tries < 5:
            timeout = 3*tries
            msg = f"Bad http response {r.status_code} from CAZy.org, retrying in {timeout} seconds..."
            logger.warning(msg)
            time.sleep(timeout)
            return self.get_clean_html_text(url_cazy, tries + 1)
        elif r.status_code != 200:
            msg = f"ERROR: HTTP response status code {r.status_code}"
            logger.error(msg)
            logger.error("Max tries to retrieve data from CAZy.org reached.")
            msg = f"CAZy.org might be down, or failing to fulfill requests to {url_cazy}"
            logger.error(msg)
            logger.error("Aborting current pipeline iteration...")
            raise PipelineException(f"Bad HTTP response code {r.status_code} from {url_cazy}, max tries exceeded.")

        # remove bad tag spacing and other tag fixes
        clean_text = self.fix_spaces.sub("><", r.text)
        clean_text = self.fix_white.sub("#ffffff", clean_text)

        return clean_text


def cazy_query(family, cazy_folder, scrape_mode, get_fragments, verbose, domain_mode,
               logger: Logger = getLogger("PipelineLogger")):
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

    cazymes: dict[str, CazymeMetadataRecord] = {}
    genbank_duplicates = []

    # loop through all pages of characterized CAZymes for selected family
    while counted < count:
        if counted > 0:
            # load new page on second iteration and above
            url = url_cazy + f"?debut_FUNC={str(counted)}#pagination_FUNC"
            clean_text = html_get.get_clean_html_text(url, logger=logger)
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

            if valid_genbank_gene(genbank, verbose) and genbank not in cazymes and \
                    (get_fragments or not protein_name.__contains__("fragment")):
                # todo: change cazymes from a dict of lists to a dict of dicts (or a dict of Namespace objects? or dict
                #  of custom class?) to get named json categories in output. THIS WILL BREAK DATA IMPORT INTO R SCRIPT
                # cazymes[genbank] = [protein_name, ec_num, org_name, None, uniprot, pdb, subfamily]  # None is domain, filled later
                cazymes[genbank] = CazymeMetadataRecord(protein_name=protein_name, ec_num=ec_num, org_name=org_name,
                                                        uniprot=uniprot, pdb=pdb, family=family,
                                                        classfamily=family.split('_')[0], subfamily=subfamily,
                                                        genbank=genbank, protein_id=genbank)
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
                    try:
                        classfam, subfam = cazyme_class.split('_')
                    except ValueError:
                        classfam = cazyme_class
                        subfam = None
                    # genbank_query = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/" + )

                    # all_cazymes[genbank] = [f"Uncharacterized {cazyme_class}", None, organism, domain, None, None, None]
                    all_cazymes[genbank] = CazymeMetadataRecord(protein_name=f"Uncharacterized {cazyme_class}",
                                                                org_name=organism, domain=domain, protein_id=genbank,
                                                                genbank=genbank, classfamily=classfam, subfamily=subfam,
                                                                family=cazyme_class)
                # we check genbank not in cazymes to prevent reporting characterized as duplicates
                elif genbank not in cazymes:
                    uncharacterized_duplicate += 1
                    msg = f"""DUPLICATE - UNCHARACTERIZED CAZYME NOT ADDED:\n
                              Organism: {organism}\n
                              Genbank: {genbank}\n
                              File Line #: {entry_reader.line_num}\n
                              \n"""
                    logger.debug(msg)
                    if verbose:
                        print(msg)
                elif genbank in cazymes:
                    # add domain to characterized entry
                    unchar_char_duplicate += 1
                    cazymes[genbank].domain = domain
                else:
                    msg = f"""UNCHARACTERIZED CAZYME NOT PARSED CORRECTLY:\n
                            Row data: {row}
                            Organism: {organism}\n
                            Genbank: {genbank}\n
                            File Line #: {entry_reader.line_num}\n
                            \n"""
                    logger.error(msg)
                    logger.error("Uncharacterized entry not parsed correctly, please report this as a bug on the "
                                 "SACCHARIS github issue tracker.")
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
                    cazymes[genbank].domain = domain

    # Filter for correct domain
    wrong_domain_characterized = 0
    wrong_domain_uncharacterized = 0
    domain_cazymes: dict[str, CazymeMetadataRecord] = {}
    for item in cazymes:
        # bitwise comparison against bitmask
        if cazymes[item].domain and Domain[cazymes[item].domain.upper()].value & domain_mode:
            domain_cazymes[item] = cazymes[item]
        else:  # invalid domain
            if cazymes[item].ec_num or cazymes[item].uniprot or cazymes[item].pdb:
                wrong_domain_characterized += 1
            else:
                wrong_domain_uncharacterized += 1

    # # This list comprehension works fine, but I changed to the for loop above to track wrong_domain vars
    # domain_cazymes = {item for item in cazymes if Domain[cazymes[item][3].upper()].value & domain_mode}
    # wrong_domain = len(cazymes) - len(domain_cazymes)

    cazy_added -= wrong_domain_characterized
    uncharacterized_added -= wrong_domain_uncharacterized

    logger.debug("Done retrieving cazymes data from CAZy!")
    print("Done retrieving cazymes data from CAZy!")
    # stats = [cazy_retrieved, cazy_added, cazy_duplicate, cazy_fragments, cazy_missing,
    #          uncharacterized_retrieved, uncharacterized_added, uncharacterized_duplicate]
    stats = [cazy_retrieved, cazy_added, cazy_duplicate, cazy_fragments, cazy_missing, wrong_domain_characterized,
             uncharacterized_retrieved, uncharacterized_added, uncharacterized_duplicate, wrong_domain_uncharacterized]

    if stats[0] != stats[1] + stats[2] + stats[3] + stats[4] + stats[5] or stats[6] != stats[7] + stats[8] + stats[9]:
        logger.warning("Summary statistics on CAZy retrieval do not add up, the statistics are not reliable, "
                       "please file a bug report with the developer.")

    return domain_cazymes, stats


def main(family, output_folder: str | os.PathLike, scrape_mode, get_fragments=False, verbose=False, force_update=False,
         ncbi_query_size=200, domain_mode=0b11111, skip_ask=False, logger: Logger = getLogger()):
    api_key, ncbi_email, ncbi_tool = load_from_env(skip_ask=skip_ask)
    # Folder and file output setup
    fasta_file = os.path.join(output_folder, f"{family}_{scrape_mode.name}_cazy.fasta")
    data_file = os.path.join(output_folder, f"{family}_{scrape_mode.name}_data.json")
    stats_file = os.path.join(output_folder, f"{family}_{scrape_mode.name}_stats.json")

    try:
        if os.path.isfile(fasta_file) and not force_update:
            print("Loading data from previous run. \nIf you wish to run a new query to the CAZy database, run "
                  "SACCHARIS 2 with --fresh")
            logger.info(f"Loading data from previous run. Data file: {data_file} ; Stats file: {stats_file}")
            with open(data_file, 'r', encoding='utf-8') as f:
                cazyme_dict = json.loads(f.read())
                cazymes = {id: CazymeMetadataRecord(**record) for id, record in cazyme_dict.items()}
            with open(stats_file, 'r', encoding='utf-8') as f:
                stats = json.loads(f.read())
            return fasta_file, cazymes, stats
    except (IOError, JSONDecodeError) as e:
        logger.debug(e)
        logger.warning(f"Error reading data from previous run... Data file: {data_file} Stats file: {stats_file}")
        if skip_ask:
            logger.warning("Continuing with fresh data download...")
        else:
            answer = ask_yes_no("An error occured while loading previous data. You can continue with a fresh download "
                                "of data from CAZy and NCBI, but this will overwrite the old, possible corrupted "
                                "files. Would you like to continue with fresh data?",
                                "Continuing with fresh data...", "Exiting...")
            if answer:
                logger.debug("User chose to continue with fresh data.")
            else:
                msg = "User chose to exit pipeline due to corrupted data load from previous run."
                logger.debug(msg)
                raise PipelineException(msg)

    cazymes, cazy_stats = cazy_query(family, output_folder, scrape_mode, get_fragments, verbose, domain_mode)

    # Take the accession numbers from the dict, convert to list, and query genbank in batches
    accession_list = list(cazymes.keys())

    fasta_data, queried, retrieved = ncbi_protein_query(accession_list, api_key, ncbi_email, ncbi_tool, False, logger,
                                                        ncbi_query_size)

    seq_list = list(SeqIO.parse(StringIO(fasta_data), "fasta"))

    cazy_stats.append(queried)
    cazy_stats.append(retrieved)

    try:
        # write standard fasta file
        with open(fasta_file, 'w') as f:
            f.write(fasta_data)
        # write data file of all the ancillary data as a dict
        with open(data_file, 'w', encoding='utf-8') as f:
            json.dump(cazymes, f, default=vars, ensure_ascii=False, indent=4,)
        # write stats file
        with open(stats_file, 'w', encoding='utf-8') as f:
            json.dump(cazy_stats, f, ensure_ascii=False, indent=4)
    except IOError as e:
        logger.error("IOError:", e)
        logger.error(f"Failed to write cazyme files in output folder: {output_folder}")
        raise PipelineException(f"Cannot write cazyme data to drive. Check that you have file write permissions in the "
                                f"output folder: {output_folder}") from e

    summary_msg = ""
    summary_msg += f"Characterized entries retrieved from CAZy database: {cazy_stats[0]}\n"
    summary_msg += f"Characterized entries added to dataset: {cazy_stats[1]}\n"
    if cazy_stats[2] > 0:
        summary_msg += f"Characterized entries with duplicate accessions not added: {cazy_stats[2]}\n"
    if cazy_stats[3] > 0:
        summary_msg += f"Characterized entries that are fragments not added: {cazy_stats[3]}\n"
    if cazy_stats[4] > 0:
        summary_msg += f"Characterized entries with missing accession not added: {cazy_stats[4]}\n"
    if cazy_stats[5] > 0:
        summary_msg += f"Characterized entries with wrong domain not added: {cazy_stats[5]}\n"
    if scrape_mode == Mode.ALL_CAZYMES:
        summary_msg += f"Uncharacterized entries retrieved from CAZy database: {cazy_stats[6]}\n"
        summary_msg += f"Uncharacterized entries added to dataset: {cazy_stats[7]}\n"
        if cazy_stats[7] > 0:
            summary_msg += f"Uncharacterized entries with duplicate accessions not added: {cazy_stats[8]}\n"
        if cazy_stats[8] > 0:
            summary_msg += f"Uncharacterized entries with wrong domain not added: {cazy_stats[9]}\n"
    summary_msg += f"Total number of queried accessions: {cazy_stats[10]}\n"
    summary_msg += f"Total number of results from NCBI: {cazy_stats[11]}\n"
    print(summary_msg)
    logger.info(summary_msg)
    if cazy_stats[10] != cazy_stats[11]:
        logger.warning("MISSING DATA FROM NCBI. SOME FASTA RESULTS COULD NOT BE QUERIED!!!")
    if cazy_stats[2] > 0 or cazy_stats[3] > 0 or cazy_stats[4] > 0 or cazy_stats[7] > 0 \
            or cazy_stats[8] != cazy_stats[9]:
        print("Details on duplicate/fragment/missing accessions printed above." if verbose else
              " Enable verbose argument for more details on duplicate/fragment/missing accessions.")
    print("\n")

    return fasta_file, cazymes, cazy_stats, seq_list


if __name__ == "__main__":
    # NOTE: this section is only used for testing and debugging purposes, this is not meant to be used!!!
    print("This file is ONLY meant to be run for debugging and development purposes!!!")
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
