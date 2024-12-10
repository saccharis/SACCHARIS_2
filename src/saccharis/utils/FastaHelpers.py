import datetime
import logging
import os
import re
from io import TextIOBase
from logging import Logger

from Bio.SeqIO import parse, write
from Bio.SeqRecord import SeqRecord

from saccharis.utils.Formatting import CazymeMetadataRecord


# todo: concatenate_multiple_fasta IS FROM DEV21! ***Delete*** concatenate_multiple_fasta() once
#  parse_multiple_fasta() is fixed and working!
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


def parse_multiple_fasta(fasta_handles: list[str | os.PathLike | TextIOBase], output_folder: str | os.PathLike = None,
                         logger: Logger = logging.getLogger(), source_override: str = None) \
        -> (dict[str:SeqRecord], dict[str:CazymeMetadataRecord], str):
    """
    Takes multiple fasta files as input, concatenates them together, optionally writes data to disk as a single FASTA file, and
    returns the output file path, a dict or sequence data records and a dict of sequence metadata records. Duplicate
    record IDs across multiple files are handled by appending duplicate IDs with '[Duplicate-User-ID-*]', with * being
    the numbered duplicate (e.g. if the ID 123456789 was duplicated, the second occurrence would be appended like so
    '123456789[Duplicate-User-ID-2]'.

    :param fasta_handles: A list of FASTA files which will be read in and concatenated together.
    :param output_folder: A folder to write the output FASTA file to.
    :param logger: Optional logging object to save logs for debugging purposes. Defaults to root logger.
    :return: Output file path, dict of metadata records, dict of sequence records
    :param source_override: A string which will be set as the source_file variable in CazymeMetadatRecord objects
    instead of the path in fasta_handles the data was read from. Used internally when parsing fasta files that
    actually come from a web API (e.g. "NCBI Gene" is put in source_file when the record was downloaded using the NCBI
    datasets API).
    """
    metadata_dict: dict[str:CazymeMetadataRecord] = {}
    all_seqs: dict[str:SeqRecord] = {}
    duplicate_counts: dict[str:int] = {}

    for path in fasta_handles:
        try:
            seqs = list(parse(path, 'fasta'))
        except FileNotFoundError as err:
            msg = f"ERROR: File path \"{err.filename}\" for provided user sequences is invalid! Did you "
            f"type it correctly?"
            logger.exception(msg)
            raise UserWarning(msg) from err
        except Exception as err:
            try:
                seqs = list(parse(path, 'fasta-2line'))
            except Exception as other:
                msg = "WARNING: Unknown error occurred while parsing user sequences. User sequences not "
                "included in analysis!\nPlease check that the file format is valid."
                logger.error("Exception 1:", err.args[0])
                logger.error("Exception 2:", other.args[0])
                logger.exception(msg)
                raise UserWarning(msg) from other

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

            if len(fasta_handles) > 1 and not source_override:
                record.description += f" SACCHARIS merged record from {path}"
            species_match = re.search(r'\[(.+)\]', record.description)
            new_record = CazymeMetadataRecord(source_file=source_override if source_override else str(path),
                                              protein_id=record_id,
                                              protein_name=record.description,
                                              org_name=species_match.group(1) if species_match else None)
            metadata_dict[record_id] = new_record
            all_seqs[record_id] = record

    out_path = None
    if output_folder:
        if not os.path.exists(output_folder):
            os.mkdir(output_folder)
        filename = f"merged_user_fasta-{datetime.datetime.now().strftime('%d-%m-%y_%H-%M')}.fasta"
        out_path = os.path.join(output_folder, filename)
        write(list(all_seqs.values()), out_path, 'fasta')

    return all_seqs, metadata_dict, out_path
