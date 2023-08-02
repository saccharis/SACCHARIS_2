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

# Dependency imports
from Bio.SeqIO import parse, write

# Internal imports
from saccharis.NCBIQueries import download_proteins_from_genomes
from saccharis.utils.PipelineErrors import FileError
from saccharis.utils.PipelineErrors import UserError
from saccharis.utils.PipelineErrors import NewUserFile
from saccharis.utils.UserInput import ask_yes_no
from saccharis.utils.UserFastaRename import rename_fasta_file
from saccharis.utils.Formatting import CazymeMetadataRecord

buf_size = 1048576  # 1 megabyte chunks
dict_filename = "user_runs.json"

# user_file = re.sub('\\', "\\\\", user_file)


class UIDValidator:
    uid_checker = re.compile('U[0-9]{1,9}')

    def __init__(self):
        self.uid_list = {}
        pass

    def validate_uid(self, record_to_test):
        match = self.uid_checker.match(record_to_test.id)
        if not match:
            raise UserError(f"Record with id {record_to_test.id} does not have a valid user header.")
        if record_to_test.id in self.uid_list:
            raise UserError(f"User id {record_to_test.id} is duplicated in user file!")
        self.uid_list[record_to_test.id] = True
        return True


def delete_files(folder):
    print(f"Deleting previous user files in {folder}...")
    for root, dirs, files in os.walk(folder):
        for file in files:
            os.remove(os.path.join(root, file))


# NOTE: This function does not guarantee that the file contains valid headers! Make sure to validate header ID before
# calling!
def calculate_user_run_id(input_file, output_folder):
    #   Calculate md5 hash of the user file, so we can disambiguate multiple user file runs
    md5 = hashlib.md5()
    try:
        with open(input_file, 'rb') as f:
            data = f.read(buf_size)
            while data:
                md5.update(data)
                data = f.read(buf_size)
        print(f"MD5 of user file is: {md5.hexdigest()}")
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
            print("WARNING:", e.args[0])
            print("WARNING: Previous user run with user file detected, but data(corrupted?) was not loaded properly.\n"
                  "WARNING: Script will continue with fresh user appended files, old files are being deleted.")
            #  json corrupt? need to delete all files in user directory, since we can't match hashes to user files
            delete_files(output_folder)
            user_dict = {}
    else:
        user_dict = {}
    user_run += len(user_dict)
    if md5.hexdigest() in user_dict:
        user_run = user_dict[md5.hexdigest()]
    else:
        user_dict[md5.hexdigest()] = user_run

    try:
        with open(user_dict_path, 'w', encoding="utf-8") as f:
            json.dump(user_dict, f, ensure_ascii=False, indent=4)
    except IOError as error:
        raise FileError("Cannot write user file hash information to file.") from error

    return user_run


def concatenate_multiple_fasta(fasta_filenames: list[str], output_folder: str, logger: logging.Logger = None)\
                                                                                                        -> [str, dict, dict]:
    metadata_dict = {}
    all_seqs = []
    # if len(fasta_filenames) == 1:
    #     seqs = parse(fasta_filenames[0], 'fasta')
    #     return fasta_filenames[0], {record.id: CazymeMetadataRecord(source_file=fasta_filenames[0],
    #                                                                 protein_id=record.id,
    #                                                                 protein_name=record.name)
    #                                                                 for record in seqs}, seqs

    for file in fasta_filenames:
        seqs = parse(file, 'fasta')
        for record in seqs:
            if record.id in metadata_dict:
                raise UserError(f"Multiple input files contain record id: '{record.id}'. Please rename record ids in "
                                f"FASTA headers for uniqueness.")
            if len(fasta_filenames) == 1:
                record.description += f" SACCHARIS merged record from {file}"
            species_match = re.search(r'\[(.+)\]', record.description)
            new_record = CazymeMetadataRecord(source_file=file,
                                              protein_id=record.id,
                                              protein_name=record.description,
                                              org_name=species_match.group(1) if species_match else None)
            metadata_dict[record.id] = new_record
            all_seqs.append(record)

    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    filename = f"merged_user_fasta-{datetime.datetime.now().strftime('%d-%m-%y_%H-%M')}.fasta"
    out_path = os.path.join(output_folder, filename)
    write(all_seqs, out_path, 'fasta')

    return out_path, metadata_dict, all_seqs


def run(user_file_path, file_to_append, output_folder, verbose=False, force_update=False, auto_rename=False,
        _run_id=None):
    # Delete all files in user folder if doing a fresh run, since they may contain outdated data from CAZy/NCBI/User
    if force_update:
        delete_files(output_folder)

    # Load sequences from the user file
    try:
        parsed_user = parse(user_file_path, 'fasta')
    except FileNotFoundError as e:
        raise UserWarning(f"ERROR: File path \"{e.filename}\" for provided user sequences is invalid! Did you type "
                          f"it correctly?")
    except Exception as e:
        try:
            parsed_user = parse(user_file_path, 'fasta-2line')
        except Exception as other:
            print("Exception 1:", e.args[0])
            print("Exception 2:", other.args[0])
            # todo: add error logging?
            raise UserWarning("WARNING: Unknown error occurred while parsing user sequences. User sequences not "
                              "included in analysis!\nPlease check that the file format is valid.")

    # Validate user seqs, including checking for valid user IDs, duplicate user ids, and querying user about renaming
    # logic when the file does not meet requirements
    validator = UIDValidator()
    user_seqs = []
    try:
        for user_record in parsed_user:
            validator.validate_uid(user_record)
            user_seqs.append(user_record)
            if verbose:
                print(f"Valid User ID: {user_record.id}")
    except UserError as error:
        rename = False
        if not auto_rename:
            print(f"ERROR: {error.msg}")
            rename = ask_yes_no("Would you like to prepend your user sequences with the proper header?",
                                "Prepending user headers in new file...",
                                "Cancelling SACCHARIS pipeline...")
        if rename or auto_rename:
            new_user_path = rename_fasta_file(user_file_path)
            raise NewUserFile(new_user_path)
        else:
            sys.exit(3)

    # Note: you should not be passing in arguments to _run_id from calls outside the pipeline, since it should be filled
    # with a value based on md5 hash information
    if not _run_id:
        _run_id = calculate_user_run_id(user_file_path, output_folder)

    output_filename = re.sub(r"\.fasta", f"_UserFile{_run_id:05d}.fasta", os.path.split(file_to_append)[1])
    output_file_path = os.path.join(output_folder, output_filename)

    #  Check for zero seq user file and preexisting output
    user_seq_count = len(user_seqs)
    if user_seq_count == 0:
        raise UserWarning("ERROR: User sequences file does not contain any valid entries! \n"
                          "Please check that file format is valid!")
    if os.path.isfile(output_file_path):
        return output_file_path, user_seq_count, _run_id

    #   load family seqs
    try:
        family_seqs = list(parse(file_to_append, 'fasta'))
    except FileNotFoundError:
        raise UserWarning("ERROR: Filename for file to append is invalid!")

    #  Append user seqs to family seqs
    new_seqs = family_seqs + user_seqs

    with open(output_file_path, 'w', newline="\n") as f:
        write(new_seqs, f, 'fasta')

    return output_file_path, user_seq_count, _run_id


def cli_main():
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
    run(args.userfile, args.familyfile, args.directory, args.verbose, args.fresh)


if __name__ == '__main__':
    cli_main()
