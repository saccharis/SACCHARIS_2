###############################################################################
#
#
# SACCHARIS 2.0 author: Alexander Fraser, https://github.com/AlexSCFraser
# License: GPL v3
###############################################################################
import math
from dataclasses import dataclass
from saccharis.utils.PipelineErrors import PipelineException


@dataclass
class CazymeMetadataRecord:
    """Class for keeping track of metadata about sequences being analyzed. Metadata gets serialized into the final
    JSON file which is an output of the main pipeline and used as input to phylogeny rendering."""
    protein_id: str
    genbank: str = None
    protein_name: str = None
    org_name: str = None
    domain: str = None
    family: str = None
    subfamily: str = None
    ec_num: str = None
    uniprot: str = None
    pdb: str = None
    module_start: str = None
    module_end: str = None
    diamond_prediction: str = None
    ecami_prediction: str = None
    source_file: str = None


def make_metadata_dict(cazy_accession_dict: dict[CazymeMetadataRecord], cazyme_module_list, bounds_dict,
                       merged_dict: dict[CazymeMetadataRecord], ecami_dict: dict, diamond_dict: dict):
    new_cazyme_dict = {}
    for module in cazyme_module_list:
        if module.__contains__("<"):
            module_id = module.split("<")[0]
        else:
            module_id = module

        # todo: delete this whole if/else section, it's been replaced by loading CazymeMetadataRecord dicts and
        #  updating records
        if module_id in cazy_accession_dict:  # module came from a CAZyme downloaded from CAZy
            genbank = module_id
            try:
                protein_name, ec_num, org_name, domain, uniprot, pdb, subfamily = cazy_accession_dict[genbank]
            except ValueError:
                protein_name, ec_num, org_name, domain, uniprot, pdb = cazy_accession_dict[genbank]
                subfamily = None  # todo: remove this branch later, it's just backwards compatibility for alpha output
            source_file = "cazy_website"
        else:  # module came from a user sequence/genome
            genbank = None
            protein_name, ec_num, org_name, domain, uniprot, pdb, subfamily = \
                None, None, None, None, None, None, None
            source_file = merged_dict[module_id].source_file if module_id in merged_dict else None

        #     Note: merged dict should not typically be none. If source_file is None, we are either working with
        #     pre-release metadata(which should probably be discarded by the time anyone reads this) or something has
        #     gone wrong, so consider a change the above to simply at a later date:
        #   source_file = merged_dict[module_id]
            if merged_dict and module_id not in merged_dict:
                print(f"ERROR: bad loading of data from merged fasta dictionary. {module_id} not in merged_dict")

        ecami_prediction = ecami_dict[module] if module in ecami_dict else None
        diamond_prediction = diamond_dict[module] if module in diamond_dict else None

        if module_id in merged_dict:
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
