###############################################################################
#
#
# SACCHARIS 2.0 author: Alexander Fraser, https://github.com/AlexSCFraser
# License: GPL v3
###############################################################################
import argparse
import os
from copy import copy, deepcopy

from Bio import SeqIO, SeqRecord

from saccharis.utils.Formatting import CazymeMetadataRecord


def rename_fasta_file(source_file_path, output_file_path=None):
    if not os.path.isfile(source_file_path):
        raise UserWarning("Filename is not valid! Please check that you have entered the correct file location, "
                          "including spelling.")

    if output_file_path is None:
        folder_name, file_name = os.path.split(source_file_path)
        name_parts = file_name.split('.')
        if len(name_parts) > 1:
            name_parts[-2] += "_UserFormat"
        elif len(name_parts) == 1:
            name_parts[0] += "_UserFormat"
            name_parts.append("fasta")
        output_file_path = os.path.join(folder_name, '.'.join(name_parts))

    try:
        seqs = []
        for record in SeqIO.parse(source_file_path, 'fasta'):
            seqs.append(record)
    except ValueError as error:
        print("Error:", error.args[0])
        raise UserWarning("Source file has invalid data, check that it is in FASTA format!") from error
    except IOError as error:
        print("Error:", error.args[0])
        raise UserWarning("Source file read error! Check that you have filesystem permissions, a common source of this"
                          " error.")
    except Exception as error:
        print("ERROR:", error.args[0])
        raise UserWarning(f"Unknown error occurred while reading fasta file {source_file_path}")

    prepend_user_headers(seqs, inplace=True)

    try:
        SeqIO.write(seqs, output_file_path, "fasta")
        print(f"Successfuly wrote renamed fasta file to {output_file_path}")
    except IOError as error:
        print("ERROR:", error.args[0])
        raise UserWarning("Error occurred trying to write renamed fasta sequences to file. Check filesystem "
                          "permissions, a common source of this error.")

    return output_file_path


def prepend_user_headers(seqs: list[SeqRecord], metadata: dict[str, CazymeMetadataRecord] = None, inplace: bool = False):
    if len(seqs) == 0:
        raise UserWarning("File contains no valid sequences! Check that the file is in a valid FASTA format.")
    new_seqs = []
    new_metadata = {}
    for i, record in enumerate(seqs):
        new_id = f"U{i:09d}"
        if inplace:
            new_metadata[new_id] = metadata[record.id]
            new_record = record
        else:
            new_metadata[new_id] = deepcopy(metadata[record.id])
            new_record = deepcopy(record)

        new_record.description = f"{new_id} {record.description}"
        new_record.name = f"{new_id}"
        new_record.id = f"{new_id}"
        new_seqs.append(new_record)

    return new_seqs, new_metadata


def cli_main():
    parser = argparse.ArgumentParser(description="SACCHARIS utility to rename entries in FASTA files with appropriate "
                                                 "headers to use in the SACCHARIS pipeline.")
    parser.add_argument("source_file", type=str, help="Path to the FASTA formatted file you would like to "
                                                      "rename the entries of")
    parser.add_argument("--destination_file", "-d", type=str, help="Filename or file path that output fasta file will"
                        "be written to. Default is to place a file named <source_filename> in the same folder as the "
                        "source file.", default=None)
    args = parser.parse_args()

    try:
        rename_fasta_file(args.source_file, args.destination_file)
    except UserWarning as error:
        print("ERROR:", error.args[0])


if __name__ == "__main__":
    cli_main()
