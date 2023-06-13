###############################################################################
#
#
# SACCHARIS 2.0 author: Alexander Fraser, https://github.com/AlexSCFraser
# License: GPL v3
###############################################################################
# Internal imports
import argparse
import io
import json
import math
import os
import re
import shutil
import sys
import csv
import textwrap
from collections import defaultdict
from contextlib import redirect_stdout
# Dependency imports
from dbcan_cli import run_dbcan
# Internal imports
from saccharis.ExtractAndPruneCAZymes import download_database
from saccharis.utils.FamilyCategories import Matcher
from saccharis.utils.FamilyCategories import get_category_list
from saccharis.utils.FamilyCategories import get_user_categories
from saccharis.utils.PipelineErrors import UserError
from saccharis.utils.UserInput import ask_yes_no
from saccharis.utils.FamilyCategories import save_family_iterable_json
from saccharis.utils.AdvancedConfig import get_db_folder


def extract_families_hmmer(fasta_filepath, output_folder, threads, hmm_eval=1e-15, hmm_cov=0.35):

    download_database()
    dbcan_output = io.StringIO()
    print(f"Screening {fasta_filepath} for CAZyme modules with hmmer settings: evalue threshold {hmm_eval} and "
          f"coverage {hmm_cov}...")
    with redirect_stdout(dbcan_output):
        run_dbcan.run(fasta_filepath, "protein", outDir=output_folder, dbDir=get_db_folder(), hmm_cpu=threads,
                      tool_arg="hmmer", hmm_eval=hmm_eval, hmm_cov=hmm_cov)

    dbcan_out_file = os.path.join(output_folder, "hmmer.out")
    # family_dict = extract_hmmer_families(dbcan_out_file)

    #   Filter hmmer output for families
    with open(dbcan_out_file, 'r', newline='\n') as hmmer_tsv:
        entry_reader = csv.reader(hmmer_tsv, delimiter='\t')
        hmmer_list = list(entry_reader)
        matcher = Matcher()
        # filters family entries from output
        family_list = [matcher.extract_cazy_family(entry[0]) for entry in hmmer_list if entry[0] != "HMM Profile"]

        # creates a dict with counts of family groupings
        family_dict = defaultdict(int)
        for entry in family_list:
            family_dict[entry] += 1
            if entry.__contains__('_'):
                family_dict[entry[0: entry.find('_')]] += 1

    return family_dict


def get_user_selection(family_dict):
    max_fam_length = max(len(key) for key in family_dict.keys())
    max_num_length = max(len(str(num)) for num in family_dict.values())
    tab_count = math.ceil((max_fam_length + max_num_length + 1)/4)
    entry_width = tab_count * 4
    console_width = shutil.get_terminal_size()[0]
    entry_count = math.floor(console_width / entry_width)
    # families_string = "\t".join(f"{key:{max_fam_length}s}:{value:<{max_num_length}d}"
    #                             for key, value in family_dict.items())
    families_string = "".join(f"{f'{key}:{value}':{entry_width}s}" for key, value in family_dict.items())
    print("\nCounts of the following CAZyme families and/or subfamilies were found in the provided user sequences:")
    print(textwrap.fill(families_string, width=entry_count*entry_width))
    # print(family_dict.items())
    try:
        user_groups = input("Please enter a space separated list of the above families and/or SACCHARIS family "
                            "categories you would like to run the pipeline on:\n")
    except KeyboardInterrupt:
        print('\n')
        sys.exit(2)
    user_groups = user_groups.split(' ')
    fam_check = Matcher()
    user_list = []
    for item in user_groups:
        test_item = item.upper()
        if fam_check.valid_cazy_family(test_item):
            if test_item in family_dict:
                user_list.append(test_item)
            else:
                print(f"{test_item} is a valid family, but not found in the user sequences. Please remove it from your "
                      f"list and try again.")
                return get_user_selection(family_dict)
        else:
            try:
                category = get_category_list(item)
                for family in category:
                    if family in family_dict:
                        user_list.append(family)

            except UserError:
                print(f"{test_item} is neither a valid family or a family category. Check spelling and try submitting "
                      f"the list again.")
                return get_user_selection(family_dict)

    return user_list


def choose_families_from_fasta(fasta_filepath, output_folder, threads):
    # run_dbcan.run(fasta_filepath, "protein", outDir=output_folder, dbDir=folder_db, hmm_cpu=threads, tool_arg="hmmer")

    family_dict = extract_families_hmmer(fasta_filepath, output_folder, threads)

    if len(family_dict) < 1:
        raise UserError("No cazymes in selected families found in user fasta sequences!")
    family_list = get_user_selection(family_dict)
    return family_list


def cli_cazome():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_file", "-f", type=str, help="The fasta file you wish to screen a cazome of.",
                        required=True)
    parser.add_argument("--out_folder", "-o", type=str, help="The folder to output hmmer and JSON results to. "
                        "Defaults to current working directory.", default=os.getcwd())
    parser.add_argument("--threads", "-t", type=int, default=math.ceil(os.cpu_count()*0.75),
                        help="HMMER allows the use of multi-core processing.  Set a number in here from"
                             " 1 to <max_cores>. The default is set at 3/4 of the number of logical cores reported by "
                             "your operating system. This is to prevent lockups if other programs are running.",
                        choices=range(1, os.cpu_count()+1))
    parser.add_argument("--family_categories", "-c", type=str, nargs='+',
                        help="A space separated list of family_categories that are in the SACCHARIS category list, "
                             "including user defined categories. When this option is included, an additional output "
                             "file will be created showing the counts of all found cazymes for each family in "
                             "family categories given to this argument. If the word all is passed in, e.g. "
                             "\'--family_categories all\', all defined categories will be calculated and output to the "
                             "file.")
    parser.add_argument('--hmm_eval', default=1e-15, type=float, help='HMMER E Value')
    parser.add_argument('--hmm_cov', default=0.35, type=float, help='HMMER Coverage val')
    args = parser.parse_args()

    input_fasta = args.fasta_file
    out_dir = args.out_folder
    cats_to_print = args.family_categories

    family_dict = extract_families_hmmer(input_fasta, out_dir, args.threads, args.hmm_eval, args.hmm_cov)

    found_file = os.path.join(out_dir, re.sub(r"\.fa.*", "_families.json", os.path.basename(input_fasta)))

    save_family_iterable_json(family_dict, found_file)

    print(f"Wrote counts of found cazyme family modules to {found_file}")

    if cats_to_print:
        user_categories = get_user_categories()
        categories = {}

        if "all" in cats_to_print or "ALL" in cats_to_print:
            for category_name, family_list in user_categories.items():
                counts = {family: (family_dict[family] if family in family_dict else 0) for family in family_list}
                categories[category_name] = counts
        else:
            for cat in cats_to_print:
                if cat in user_categories:
                    counts = {family: (family_dict[family] if family in family_dict else 0)
                              for family in user_categories[cat]}
                    categories[cat] = counts
                else:
                    print(f"\'{cat}\' category not found in family categories. Check spelling or add this "
                          f"family category to get the formatted results.")
                    skip = ask_yes_no(f"Continue anyway and skip {cat} category?", "Continuing...", "Exiting...")
                    if not skip:
                        sys.exit()

        category_file = os.path.join(out_dir, re.sub(r"\.fa.*", "_categories.json", os.path.basename(input_fasta)))
        with open(category_file, 'w', encoding="utf-8") as jsonfile:
            json.dump(categories, jsonfile, ensure_ascii=False, indent=4)

        print(f"Wrote counts of found cazyme family modules in spcified categories to {category_file}")


if __name__ == "__main__":
    cli_cazome()
