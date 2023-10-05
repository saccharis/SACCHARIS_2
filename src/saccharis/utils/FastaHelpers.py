import datetime
import os
import re
from io import TextIOBase
from logging import Logger

from Bio.SeqIO import parse, write
from Bio.SeqRecord import SeqRecord

from saccharis.utils.Formatting import CazymeMetadataRecord
from saccharis.utils.PipelineErrors import UserError


def parse_multiple_fasta(fasta_handles: list[str | os.PathLike | TextIOBase], output_folder: str | os.PathLike = None,
                         logger: Logger = None, source_override: str = None) \
        -> (list[SeqRecord], dict[str:CazymeMetadataRecord], str):

    metadata_dict: dict[str:CazymeMetadataRecord] = {}
    all_seqs: list[SeqRecord] = []

    for path in fasta_handles:
        try:
            seqs = list(parse(path, 'fasta'))
        except FileNotFoundError as err:
            raise UserWarning(f"ERROR: File path \"{err.filename}\" for provided user sequences is invalid! Did you "
                              f"type it correctly?") from err
        except Exception as err:
            try:
                seqs = list(parse(path, 'fasta-2line'))
            except Exception as other:
                logger.error("Exception 1:", err.args[0])
                logger.error("Exception 2:", other.args[0])
                raise UserWarning("WARNING: Unknown error occurred while parsing user sequences. User sequences not "
                                  "included in analysis!\nPlease check that the file format is valid.") from other

        for record in seqs:
            if record.id in metadata_dict:
                raise UserError(f"Multiple input files contain record id: '{record.id}'. Please rename record ids in "
                                f"FASTA headers for uniqueness.")
            if len(fasta_handles) > 1 and not source_override:
                record.description += f" SACCHARIS merged record from {path}"
            species_match = re.search(r'\[(.+)\]', record.description)
            new_record = CazymeMetadataRecord(source_file=source_override if source_override else path,
                                              protein_id=record.id,
                                              protein_name=record.description,
                                              org_name=species_match.group(1) if species_match else None)
            metadata_dict[record.id] = new_record
            all_seqs.append(record)

    out_path = None
    if output_folder:
        if not os.path.exists(output_folder):
            os.mkdir(output_folder)
        filename = f"merged_user_fasta-{datetime.datetime.now().strftime('%d-%m-%y_%H-%M')}.fasta"
        out_path = os.path.join(output_folder, filename)
        write(all_seqs, out_path, 'fasta')

    return all_seqs, metadata_dict, out_path
