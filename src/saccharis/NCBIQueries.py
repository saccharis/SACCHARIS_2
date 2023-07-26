import io
import logging
import os
import re
import time
from zipfile import ZipFile

import Bio.SeqIO
from Bio.SeqRecord import SeqRecord
import requests
from bs4 import BeautifulSoup
from ncbi.datasets import GenomeApi

from saccharis.Cazy_Scrape import NCBIException
from saccharis.utils.PipelineErrors import PipelineException

NCBI_DEFAULT_DELAY = 0.3  # this is a delay time for ncbi queries. Without it, results may be incomplete


def valid_ncbi_genome(string_to_check: str, verbose: bool = False):
    # The goal is to validate that the string matches some kind of genome accession identifier from NCBI. This page
    # describes several formats they can take, but this is incompelte because it doesn't provide detailed information on
    # the RefSeq identifiers. https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/
    # RefSeq ids are supposed to be described here, but this is also incompletel because it doesn't include 3 letter
    # codes before the underscore. https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly/
    # A very common format for reference genomes is:
    #  [three-letter alphabetical prefix]_[9 digits][.][version number]
    #   Where the three letters are GCA or GCF for GenBank and RefSeq accessions respectively.
    #   e.g. GCA_018292165.1 or

    if string_to_check is None or string_to_check == "":
        return False

    if re.match(r"((GCA)|(GCF))_\w{9}\.\d+", string_to_check):
        return True

    if verbose:
        print(f"\"{string_to_check}\" is not a valid genbank gene/protein accession identifier. If this is in error "
              f"(i.e. this string is in fact a valid genbank accession for a SINGLE gene/protein), please report this "
              f"as a bug to the developer via the github page and/or email.")
    return False


def valid_genbank_gene(string_to_check: str, verbose: bool = False):
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
        print(f"\"{string_to_check}\" is not a valid genbank gene/protein accession identifier. If this is in error "
              f"(i.e. this string is in fact a valid genbank accession for a SINGLE gene/protein), please report this "
              f"as a bug to the developer via the github page and/or email.")
    return False


# todo: call this function from gui, lol
def download_proteins_from_genomes(genome_list: list[str], out_dir: str = None, logger: logging.Logger = None,
                                   fresh: bool = False) -> (list[SeqRecord], dict[str:str]):
    seqs = []
    source_dict = {}
    api = GenomeApi()
    if out_dir:
        if not fresh:
            # todo: check for local seqs to load from each genome  instead of downloading from NCBI, if fresh == false
            pass
        outpath = os.path.join(out_dir, "ncbi_dataset.zip")
        handle = api.download_assembly_package(genome_list, include_annotation_type=["PROT_FASTA"], filename=outpath)
    else:
        handle = api.download_assembly_package(genome_list, include_annotation_type=["PROT_FASTA"])
    try:
        with ZipFile(handle) as myzip:
            for genome_id in genome_list:
                with io.TextIOWrapper(myzip.open(f"ncbi_dataset/data/{genome_id}/protein.faa"), encoding="utf-8") as myfile:
                    genome_seqs = Bio.SeqIO.parse(myfile, "fasta")
                    source_dict += dict.fromkeys(map(lambda seq: seq.id, genome_seqs), f"NCBI Genome:{genome_id}")
                    # todo: save seqs locally for later if out_dir is given
                    seqs += genome_seqs
    except Exception as e:
        if logger:
            msg = "Problem reading genome zip file downloaded from NCBI."
            logger.error(e.__traceback__)
            logger.error(msg)
            raise PipelineException(msg) from e
    finally:
        handle.close()
        os.remove(handle.name)

    return seqs, source_dict


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
    #  when the NCBI server has errors, so I'm not sure which specific exceptions to catch.
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
    # For more information on the details of the NCBI accession fasta ID format, see the following page
    # https://ncbi.github.io/cxx-toolkit/pages/ch_demo#ch_demo.id1_fetch.html_ref_fasta
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
        print(f"WARNING: Probable parsing error on accession containing a '|' character, "
              f"Please report this as a bug to the developer/maintainer through github.")

    return fasta_out, result_count


