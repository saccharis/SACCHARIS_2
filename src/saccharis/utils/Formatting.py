###############################################################################
#
#
# SACCHARIS 2.0 author: Alexander Fraser, https://github.com/AlexSCFraser
# License: GPL v3
###############################################################################
import math


def make_metadata_dict(cazy_accession_dict, cazyme_module_list, bounds_dict, merged_dict, ecami_dict, diamond_dict):
    new_cazyme_dict = {}
    for module in cazyme_module_list:
        if module.__contains__("<"):
            module_id = module.split("<")[0]
        else:
            module_id = module

        if module_id in cazy_accession_dict:  # module came from a CAZyme downloaded from CAZy
            genbank = module_id
            try:
                protein_name, ec_num, org_name, domain, uniprot, pdb, subfamily = cazy_accession_dict[genbank]
            except ValueError:
                protein_name, ec_num, org_name, domain, uniprot, pdb = cazy_accession_dict[genbank]
                subfamily = None  # todo: remove this branch later, it's just backwards compatibility for alpha output
            try:
                ecami_prediction = ecami_dict[module]
            except KeyError:
                ecami_prediction = None
            try:
                diamond_prediction = diamond_dict[module]
            except KeyError:
                diamond_prediction = None
            source_file = "cazy_website"
        else:  # module came from a user sequence/genome
            genbank = None
            protein_name, ec_num, org_name, domain, uniprot, pdb, subfamily, ecami_prediction, diamond_prediction = \
                None, None, None, None, None, None, None, None, None
            source_file = merged_dict[module_id] if module_id in merged_dict else None

        #     Note: merged dict should not typically be none. If source_file is None, we are either working with
        #     pre-release metadata(which should probably be discarded by the time anyone reads this) or something has
        #     gone wrong, so consider change the above to simply at a later date:
        # source_file = merged_dict[module_id]
            if merged_dict and module_id not in merged_dict:
                print("ERROR: bad loading of data from merged fasta dictionary.")

        entry_dict = {"genbank": genbank, "protein_name": protein_name, "ec_num": ec_num, "org_name": org_name,
                      "domain": domain, "uniprot": uniprot, "pdb": pdb, "module_start": bounds_dict[module][0],
                      "module_end": bounds_dict[module][1], "subfamily": subfamily,
                      "diamond_prediction": diamond_prediction, "ecami_prediction": ecami_prediction,
                      "source_file": source_file}

        new_cazyme_dict[module] = entry_dict

    return new_cazyme_dict


def format_time(seconds):
    if seconds > 3600:
        return f"*\t {math.floor(seconds/3600):.0f} hours, {math.floor(seconds%3600/60):.0f} minutes, " \
               f"{seconds % 60:.0f} seconds to run"
    if seconds > 60:
        return f"*\t {math.floor(seconds/60):.0f} minutes, {seconds%60:.0f} seconds to run"

    return f"*\t {seconds:.1f} seconds to run"
