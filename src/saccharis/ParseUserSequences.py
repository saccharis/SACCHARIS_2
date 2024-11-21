###############################################################################
# Functions to parse user FASTA files for use in the SACCHARIS 2.0 pipeline.
# Author: Alexander Fraser, https://github.com/AlexSCFraser
# License: GPL v3
###############################################################################
# Built in libraries
import argparse
import json
import logging
import os
import re
import hashlib
import sys
import datetime
from collections import Counter
from io import TextIOBase, StringIO
from logging import Logger, getLogger
from typing import Iterable

# Dependency imports
from Bio.SeqIO import parse, write, SeqRecord

# Internal imports
from saccharis.NCBIQueries import download_proteins_from_genomes, download_from_genes
from saccharis.utils.AdvancedConfig import get_ncbi_folder
from saccharis.utils.FastaHelpers import parse_multiple_fasta
from saccharis.utils.PipelineErrors import FileError, PipelineException, NewUserFile
from saccharis.utils.PipelineErrors import UserError
from saccharis.utils.UserInput import ask_yes_no
from saccharis.utils.UserFastaRename import prepend_user_headers, rename_fasta_file
from saccharis.utils.Formatting import CazymeMetadataRecord, seqs_to_string


buf_size = 1048576  # 1 megabyte chunks
dict_filename = "user_runs.json"

# user_file = re.sub('\\', "\\\\", user_file)


class UIDValidator:
    uid_checker = re.compile('U[0-9]{1,9}')

    def __init__(self):
        self.uid_list = {}
        pass

    def validate_uid(self, record_to_test: SeqRecord):
        match = self.uid_checker.match(record_to_test.id)
        if not match:
            raise UserError(f"Record with id {record_to_test.id} does not have a valid user header.")
        if record_to_test.id in self.uid_list:
            raise UserError(f"User id {record_to_test.id} is duplicated in user file!")
        self.uid_list[record_to_test.id] = True
        return True


def find_duplicates(elements: Iterable):
    count = Counter(elements)
    duplicates = [item for item, freq in count.items() if freq > 1]
    return duplicates


def delete_files(folder):
    print(f"Deleting previous user files in {folder}...")
    for root, dirs, files in os.walk(folder):
        for file in files:
            os.remove(os.path.join(root, file))


def calculate_user_run_id(input_handle: str | os.PathLike | TextIOBase, output_folder, logger: logging.Logger = logging.getLogger()):
    """
    Calculates md5 hash of the user file and creates a user run ID number, so we can disambiguate multiple user file
    runs to restore partially completed runs without redoing extra computation.

    NOTE: This function DOES NOT guarantee that the file contains valid headers! Make sure to validate header ID
    before calling!

    :param input_file: The file to compute a hash from.
    :param output_folder: The folder that stores a file which contains the mapping from hash to user run ID.
    :return: Returns the user run ID number.
    """
    md5 = hashlib.md5()

    try:
        if type(input_handle) == str or type(input_handle) == os.PathLike:
            with open(input_handle, 'rb') as f:
                data = f.read(buf_size)
                while data:
                    md5.update(data)
                    data = f.read(buf_size)
        else:  # type(input_handle) == TextIOBase:
            data = input_handle.read(buf_size)
            while data:
                md5.update(data.encode())
                data = input_handle.read(buf_size)
        md5_string = md5.hexdigest()
        print(f"MD5 of user file is: {md5_string}")
    except FileNotFoundError as e:
        raise UserWarning(f"ERROR: File path \"{e.filename}\" for provided user sequences is invalid! Did you type "
                          f"it correctly?")
    except Exception as e:
        print("ERROR:", e.args[0])
        raise UserWarning("Unknown file I/O error has occurred.")

    #   load dict of previous user files to determine output filename based on md5 hash
    user_run = 0
    user_dict_path = os.path.join(output_folder, dict_filename)
    if os.path.isfile(user_dict_path):
        try:
            with open(user_dict_path, 'r', encoding='utf-8') as f:
                user_dict = json.loads(f.read())
        except Exception as e:
            logger.exception(f"Uncaught exception caused failure to load user run hash mapping from JSON file at "
                             f"{user_dict_path}.")
            logger.warning("Previous user run with user file detected, but data(corrupted?) was not loaded properly.")
            logger.warning("Script will continue with fresh user appended files, old files are being deleted.")
            #  json corrupt? need to delete all files in user directory, since we can't match hashes to user files
            # todo: this would be better if it caught specific exceptions in case we need to debug this failure mode
            delete_files(output_folder)
            user_dict = {}
    else:
        user_dict = {}
    user_run += len(user_dict)
    if md5_string in user_dict:
        user_run = user_dict[md5_string]
    else:
        user_dict[md5_string] = user_run

    try:
        with open(user_dict_path, 'w', encoding="utf-8") as f:
            json.dump(user_dict, f, ensure_ascii=False, indent=4)
    except IOError as error:
        msg = f"Cannot write user file hash information to file: {user_dict_path}"
        logger.exception(msg)
        raise FileError(msg) from error

    return user_run, md5_string


# def parse_user_files(user_file_paths: list[str | os.PathLike], logger: Logger = getLogger()) \
#         -> (list[SeqRecord], dict[str:str]):
#     seqs: list[SeqRecord] = []
#     source_file: dict[str] = {}
#
#     for path in user_file_paths:
#         try:
#             parsed_user = list(parse(path, 'fasta'))
#         except FileNotFoundError as e:
#             raise UserWarning(f"ERROR: File path \"{e.filename}\" for provided user sequences is invalid! Did you "
#                               f"type it correctly?") from e
#         except Exception as e:
#             try:
#                 parsed_user = list(parse(path, 'fasta-2line'))
#             except Exception as other:
#                 logger.error("Exception 1:", e.args[0])
#                 logger.error("Exception 2:", other.args[0])
#                 raise UserWarning("WARNING: Unknown error occurred while parsing user sequences. User sequences not "
#                                   "included in analysis!\nPlease check that the file format is valid.") from other
#
#         for record in parsed_user:
#             if record.id in source_file:
#                 raise UserWarning(f"Sequence identifier {record.id} is present multiple times in your FASTA files! "
#                                   f"Please use unique sequence identifiers within and across your FASTA files.")
#             seqs.append(record)
#             source_file[record.id] = path
#
#     return seqs, source_file

# todo: concatenate_multiple_fasta IS FROM DEV21, AND WAS ADJUSTED TO FIX DUPLICATE IDS ACROSS USER FILES I NEED TO
#  TAKE THE IDEA HERE AND BRING IT INTO THE REFACTORED DEV22 METHOD OF MERGING DATA SOURCES,
#  THIS CODE MIGHT BE BETTER IN FastaHelpers MODULE
def concatenate_multiple_fasta(fasta_filenames: list[str | os.PathLike], output_folder: str | os.PathLike,
                               logger: logging.Logger = logging.getLogger()) -> [str, dict, dict]:
    """
    Takes multiple fasta files as input, concatenates them together, writes data to disk as a single FASTA file, and
    returns the output file path, a dict or sequence data records and a dict of sequence metadata records. Duplicate
    record IDs across multiple files are handled by appending duplicate IDs with '[Duplicate-User-ID-*]', with * being
    the numbered duplicate (e.g. if the ID 123456789 was duplicated, the second occurrence would be appended like so
    '123456789[Duplicate-User-ID-2]'.

    :param fasta_filenames: A list of FASTA files which will be read in and concatenated together.
    :param output_folder: A folder to write the output FASTA file to.
    :param logger: Optional logging object to save logs for debugging purposes. Defaults to root logger.
    :return: Output file path, dict of metadata records, dict of sequence records
    """
    metadata_dict = {}
    all_seqs = {}
    duplicate_counts = {}
    # if len(fasta_filenames) == 1:
    #     seqs = parse(fasta_filenames[0], 'fasta')
    #     return fasta_filenames[0], {record.id: CazymeMetadataRecord(source_file=fasta_filenames[0],
    #                                                                 protein_id=record.id,
    #                                                                 protein_name=record.name)
    #                                                                 for record in seqs}, seqs

    for file in fasta_filenames:
        seqs = parse(file, 'fasta')
        for record in seqs:
            old_id = record.id
            if record.id in metadata_dict:
                # rename existing duplicated ID
                duplicate_counts[old_id] = 1
                new_id = metadata_dict[old_id].protein_id + "[Duplicate-User-ID-1]"
                metadata_dict[old_id].protein_id = new_id
                metadata_dict[new_id] = metadata_dict[old_id]
                metadata_dict.pop(old_id)

                all_seqs[old_id].id = new_id
                all_seqs[old_id].description = all_seqs[old_id].description.replace(old_id, new_id)
                all_seqs[old_id].name = all_seqs[old_id].name.replace(old_id, new_id)
                all_seqs[new_id] = all_seqs[old_id]
                all_seqs.pop(old_id)
                logger.debug(f"User sequence with duplicate ID: {record.id} renamed to {new_id}")

            if record.id in duplicate_counts:
                duplicate_counts[record.id] += 1
                new_id = record.id + f"[Duplicate-User-ID-{duplicate_counts[record.id]}]"
                record_id = new_id
                record.id = record_id
                record.name = record.name.replace(old_id, new_id)
                record.description = record.description.replace(old_id, new_id)
                logger.debug(f"User sequence with duplicate ID: {old_id} named to {new_id}")
            else:
                record_id = record.id

            if len(fasta_filenames) == 1:
                record.description += f" SACCHARIS merged record from {file}"
            species_match = re.search(r'\[(.+)\]', record.description)
            new_record = CazymeMetadataRecord(source_file=file,
                                              protein_id=record_id,
                                              protein_name=record.description,
                                              org_name=species_match.group(1) if species_match else None)
            metadata_dict[record_id] = new_record
            all_seqs[record_id] = record

    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    filename = f"merged_user_fasta-{datetime.datetime.now().strftime('%d-%m-%y_%H-%M')}.fasta"
    out_path = os.path.join(output_folder, filename)
    write(list(all_seqs.values()), out_path, 'fasta')

    return out_path, metadata_dict, all_seqs



def merge_data_sources(cazy_seqs: list[SeqRecord] | None, cazy_metadata: dict[str:CazymeMetadataRecord] | None,
                       user_file_paths: list[str | os.PathLike], ncbi_genomes: list[str] | None,
                       ncbi_genes: list[str] | None, output_folder: str | os.PathLike, output_prefix: str = None,
                       verbose: bool = False, force_update: bool = False, auto_rename: bool = False,
                       logger: Logger = getLogger(), skip_user_ask: bool = False, ask_func=ask_yes_no):
    # Delete all files in output folder if doing a fresh run, since they may contain outdated data from CAZy/NCBI/User
    if force_update:
        delete_files(output_folder)

    if cazy_seqs is None:
        cazy_seqs = []
    if cazy_metadata is None:
        cazy_metadata = []

    non_cazy_seqs: list[SeqRecord] = []
    non_cazy_metadata: dict[str:CazymeMetadataRecord] = {}

    if user_file_paths is not None:
        try:
            fasta_seqs, fasta_metadata, _ = parse_multiple_fasta(user_file_paths, logger=logger)
            if len(fasta_seqs) == 0:
                raise UserWarning("ERROR: User sequence file(s) do not contain any valid entries! \n"
                                  "Please check that file format is valid!")
            non_cazy_seqs += fasta_seqs
            non_cazy_metadata |= fasta_metadata
        except (UserWarning, UserError) as error:
            logger.warning(error.args[0])
            answer = ask_func("Would you like to continue anyway, without user sequences in the analysis?",
                              "Continuing without user fasta sequences...",
                              "Cancelling SACCHARIS pipeline...")
            if not skip_user_ask and answer:
                print("Continuing...")
                logger.info("Continuing...")
            else:
                print("Exiting...")
                logger.info("Exiting...")
                sys.exit()

    if ncbi_genomes is not None:
        genome_seqs, genome_metadata = download_proteins_from_genomes(ncbi_genomes, out_dir=get_ncbi_folder(),
                                                                      logger=logger)
        non_cazy_seqs += genome_seqs
        non_cazy_metadata |= genome_metadata

    if ncbi_genes is not None:
        gene_seqs, gene_metadata = download_from_genes(ncbi_genes, seq_type="protein", out_dir=get_ncbi_folder(),
                                                       logger=logger, fresh=force_update)
        non_cazy_seqs += gene_seqs
        non_cazy_metadata |= gene_metadata

    # Validate user seqs, including checking for valid user IDs, duplicate user ids, and querying user about renaming
    # logic when the file does not meet requirements
    validator = UIDValidator()
    user_seqs = []
    try:
        for user_record in non_cazy_seqs:
            validator.validate_uid(user_record)
            user_seqs.append(user_record)
            if verbose:
                print(f"Valid User ID: {user_record.id}")
    except UserError as error:
        rename = False
        if not auto_rename and not skip_user_ask:
            logger.error(f"{error.msg}")
            rename = ask_func("Would you like to prepend your user sequences with the proper header?",
                              "Prepending user headers in new file...",
                              "Cancelling SACCHARIS pipeline...")
        if rename or auto_rename or skip_user_ask:
            prepend_user_headers(non_cazy_seqs, non_cazy_metadata, inplace=True)
        else:
            sys.exit(3)

    all_metadata: dict[str:CazymeMetadataRecord] = cazy_metadata | non_cazy_metadata
    all_seqs: list[SeqRecord] = cazy_seqs + non_cazy_seqs

    if len(all_seqs) != len(all_metadata):
        duplicate_ids = find_duplicates(map(lambda seq: seq.id, all_seqs))
        raise PipelineException("Length of user sequence array and user sequence metadata array do not match. This "
                                "error occurs because there are duplicate accession IDs across user files, CAZy "
                                f"sequences, and NCBI sequences. Duplicated IDs: {duplicate_ids}")

    all_seq_data = StringIO(seqs_to_string(all_seqs))
    _run_id, md5_hash = calculate_user_run_id(all_seq_data, output_folder)
    if output_prefix is None:
        output_filename = f"{md5_hash}.faa"
    else:
        output_filename = f"{output_prefix}{f'_UserFile{_run_id:05d}.fasta'}.faa"
    # output_filename = re.sub(r"\.fasta", f"_UserFile{_run_id:05d}.fasta", os.path.split(file_to_append)[1])
    output_file_path = os.path.join(output_folder, output_filename)

    #  Check for preexisting output
    if os.path.isfile(output_file_path):
        return all_seqs, all_metadata, output_file_path, _run_id

    with open(output_file_path, 'w', newline="\n") as f:
        write(all_seqs, f, 'fasta')

    return all_seqs, all_metadata, output_file_path, _run_id


def cli_main():
    # todo: delete this CLI interface? It it outdated but also unused and doesn't seem particularly useful
    parser = argparse.ArgumentParser()
    # user_file_path, file_to_append, output_folder, verbose=False, force_update=False
    parser.add_argument("--familyfile", "-f", type=str, help="These are the cazyme family sequences which user "
                        "sequences will be appended to. Sequences MUST be in FASTA FORMAT - if they are not the script "
                        "will fail. Make sure to include path with filename.", required=True)
    parser.add_argument("--userfile", "-u", type=str, help="These are the user sequences being appended to the CAzyme "
                        "family sequences.  Sequences MUST be in FASTA FORMAT - if they are not the script will "
                        "fail.  Make sure to include path with filename.", required=True)
    parser.add_argument('--directory', "-o", type=str, default=os.path.join(os.getcwd(), "output"), help='You can set '
                        'a predefined output directory with this flag, either a full path or a subfolder of the CWD.  '
                        'Default is <Current Working Directory (CWD)>/output. If you specify an absolute file path the '
                        'end directory will be used. If you specify a relative file path(e.g. just a folder name), it'
                        ' will be a subdirectory of the CWD.')
    parser.add_argument("--fresh", "-n", action="store_true", help="This is a boolean value flag that by default "
                        "is set to False, which means existing data will be reused to speed up analysis. When included,"
                        " this options forces analyses to be performed again. Saved data "
                        "from previous runs with the same filenames will be overwritten.")
    parser.add_argument("--verbose", "-v", action="store_true", help="This is a boolean value flag that by default "
                        "is set to False, which means verbose output is hidden. If you would like verbose output"
                        ", include this flag in your call.")
    args = parser.parse_args()

    # run(args.userfile, args.familyfile, args.directory, args.verbose, args.fresh)


# if __name__ == '__main__':
#     cli_main()
