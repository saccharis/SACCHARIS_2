###############################################################################
#
#
# SACCHARIS 2.0 author: Alexander Fraser, https://github.com/AlexSCFraser
# License: GPL v3
###############################################################################
import math
import subprocess
from dataclasses import dataclass
from logging import Logger, getLogger
from typing import Optional

from saccharis.utils.PipelineErrors import PipelineException


def convert_path_wsl(path: str):
    if not path.__contains__(' '):
        # this line feels unnecessary, but for some reason this function breaks on windows paths without spaces when
        # single quoted, and ALSO breaks on windows paths WITH spaces when single quoted, so we only single quote on
        # paths without spaces. Very weird behaviour. Double quotes don't work at all, they remove all the slahes.
        path = "'" + path + "'"
    return subprocess.run([f"wsl", "wslpath", path], capture_output=True, check=True).stdout.decode().strip()


@dataclass
class CazymeMetadataRecord:
    """Class for keeping track of metadata about sequences being analyzed. Metadata gets serialized into the final
    JSON file which is an output of the main pipeline and used as input to phylogeny rendering."""
    protein_id: str
    genbank: Optional[str] = None
    protein_name: Optional[str] = None
    org_name: Optional[str] = None
    domain: Optional[str] = None
    family: Optional[str] = None
    subfamily: Optional[str] = None
    ec_num: Optional[str] = None
    uniprot: Optional[str] = None
    pdb: Optional[str] = None
    module_start: Optional[str] = None
    module_end: Optional[str] = None
    diamond_prediction: Optional[str] = None
    ecami_prediction: Optional[str] = None
    source_file: Optional[str] = None


def rename_header_ids(new_user_fasta_file: str, metadata_dict: dict[str, CazymeMetadataRecord]) \
                                                                                    -> dict[str, CazymeMetadataRecord]:
    new_metadata_dict: dict[str, CazymeMetadataRecord] = {}
    with open(new_user_fasta_file, 'r') as new_user_file:
        headers = filter(lambda line: line.__contains__('>'), new_user_file.readlines())
        rename_dict = {line.split(' ')[1]: line.split(' ')[0][1:] for line in headers}
    for record_id in metadata_dict:
        new_metadata_dict[rename_dict[record_id]] = metadata_dict[record_id]

    return new_metadata_dict


def make_metadata_dict(cazy_accession_dict: dict[CazymeMetadataRecord], cazyme_module_list, bounds_dict,
                       merged_dict: dict[CazymeMetadataRecord], ecami_dict: dict, diamond_dict: dict,
                       logger: Logger = getLogger()):
    new_cazyme_dict = {}
    for module in cazyme_module_list:
        if module.__contains__("<"):
            module_id = module.split("<")[0]
        else:
            module_id = module

        # # todo: delete this whole if/else section, it's been replaced by loading CazymeMetadataRecord dicts and
        # #  updating records
        # if module_id in cazy_accession_dict:  # module came from a CAZyme downloaded from CAZy
        #     genbank = module_id
        #     try:
        #         protein_name, ec_num, org_name, domain, uniprot, pdb, subfamily = cazy_accession_dict[genbank]
        #     except ValueError:
        #         protein_name, ec_num, org_name, domain, uniprot, pdb = cazy_accession_dict[genbank]
        #         subfamily = None  # todo: remove this branch later, it's just backwards compatibility for alpha output
        #     source_file = "cazy_website"
        # else:  # module came from a user sequence/genome
        #     genbank = None
        #     protein_name, ec_num, org_name, domain, uniprot, pdb, subfamily = \
        #         None, None, None, None, None, None, None
        #     source_file = merged_dict[module_id].source_file if module_id in merged_dict else None

        #     Note: merged dict should not typically be none. If source_file is None, we are either working with
        #     pre-release metadata(which should probably be discarded by the time anyone reads this) or something has
        #     gone wrong, so consider a change the above to simply at a later date:
        #   source_file = merged_dict[module_id]
            if merged_dict and module_id not in merged_dict and module_id not in cazy_accession_dict:
                logger.error(f"Bad loading of data from merged fasta dictionary. {module_id} not in merged_dict")

        ecami_prediction = ecami_dict[module] if module in ecami_dict else None
        diamond_prediction = diamond_dict[module] if module in diamond_dict else None

        if merged_dict and module_id in merged_dict:
            merged_dict[module_id].ecami_prediction = ecami_prediction
            merged_dict[module_id].diamond_prediction = diamond_prediction
            merged_dict[module_id].module_start = bounds_dict[module][0]
            merged_dict[module_id].module_end = bounds_dict[module][1]
            entry_item = merged_dict[module_id]
        elif module_id in cazy_accession_dict:
            cazy_accession_dict[module_id].ecami_prediction = ecami_prediction
            cazy_accession_dict[module_id].diamond_prediction = diamond_prediction
            cazy_accession_dict[module_id].module_start = bounds_dict[module][0]
            cazy_accession_dict[module_id].module_end = bounds_dict[module][1]
            entry_item = cazy_accession_dict[module_id]
        else:
            raise PipelineException(f"Error in make_meatdata_dict method, it failed to receive a CazymeMetadataRecord "
                                    f"for accession id {module_id} in it's arguments")
            # todo: delete this whole code for instantiating new CazymeMetadataRecord objects here, it shouldn't run.
            #  I am keeping the code in for now as reference until loading data from the cazy_accession_dict and
            #  merged_dict is bug free and produces correct JSON files
            entry_item = CazymeMetadataRecord(protein_name=protein_name,
                                               # todo: add protein id field to capture protein id from user seqs. will be equal to
                                               #  genbank for cazymes from cazy
                                               protein_id=module_id,
                                               genbank=genbank,
                                               org_name=org_name,
                                               domain=domain,
                                               # todo: add family field to track family. This is sort of redundnant for simple
                                               #  analyses, which is why it wasn't originally included
                                               # family=family,
                                               subfamily=subfamily,
                                               ec_num=ec_num,
                                               uniprot=uniprot,
                                               pdb=pdb,
                                               module_start=bounds_dict[module][0],
                                               module_end=bounds_dict[module][1],
                                               diamond_prediction=diamond_prediction,
                                               ecami_prediction=ecami_prediction,
                                               source_file=source_file
                                              )
        # old encoding was a dict, delete this once CazymeRecord implementation is tested and working
        # entry_dict = {"genbank": genbank, "protein_name": protein_name, "ec_num": ec_num, "org_name": org_name,
        #           "domain": domain, "uniprot": uniprot, "pdb": pdb, "module_start": bounds_dict[module][0],
        #           "module_end": bounds_dict[module][1], "subfamily": subfamily,
        #           "diamond_prediction": diamond_prediction, "ecami_prediction": ecami_prediction,
        #           "source_file": source_file}

        new_cazyme_dict[module] = entry_item

    return new_cazyme_dict


def format_time(seconds):
    if seconds > 3600:
        return f"*\t {math.floor(seconds/3600):.0f} hours, {math.floor(seconds%3600/60):.0f} minutes, " \
               f"{seconds % 60:.0f} seconds to run"
    if seconds > 60:
        return f"*\t {math.floor(seconds/60):.0f} minutes, {seconds%60:.0f} seconds to run"

    return f"*\t {seconds:.1f} seconds to run"
