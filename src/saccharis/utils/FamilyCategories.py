###############################################################################
#
#
# SACCHARIS 2.0 author: Alexander Fraser, https://github.com/AlexSCFraser
# License: GPL v3
###############################################################################
# Built in libraries
import io
import json
import os
import argparse
import re
import sys
from json import JSONDecodeError
# Internal Imports
from utils.PipelineErrors import PipelineException
from utils.UserInput import ask_yes_no
from saccharis.utils.PipelineErrors import UserError
from saccharis.utils.AdvancedConfig import get_config_folder

folder_config = get_config_folder()
all_families_filename = "all_families.json"
fam_lists_filename = "family_categories.json"
default_fam_lists_file_path = os.path.join(folder_config, fam_lists_filename)

_DELETEDFAMILYLIST = ["CBM33", "CBM7", "CE10", "GH145", "GH155", "GH21", "GH40", "GH41", "GH60", "GH61", "GH69", "GT36",
                      "PL19"]


class Matcher:
    cazy_fam_regex = re.compile(r"((GH)|(PL)|(GT)|(CE)|(AA)|(CBM))\d+(_\d+)?")

    def __init__(self):
        pass

    def valid_cazy_family(self, family_string_to_test):
        return bool(self.cazy_fam_regex.fullmatch(family_string_to_test))

    def extract_cazy_family(self, string_to_extract_from):
        result = self.cazy_fam_regex.match(string_to_extract_from)
        if result:
            return result.group()
        else:
            return ""


def check_deleted_families(family):
    if family in _DELETEDFAMILYLIST:
        raise PipelineException(f"Family {family} is a deleted family. "
                                f"Check https://www.cazypedia.org/index.php/{family} for more details")


def get_default_family_categories():

    fam_lists = {"all_families": []}
    families = [("GH", 173), ("GT", 115), ("PL", 20), ("CE", 20), ("AA", 17), ("CBM", 91)]

    for fam_tuple in families:
        for num in range(1, fam_tuple[1]+1):
            if f"{fam_tuple[0]}{num}" not in _DELETEDFAMILYLIST:
                fam_lists["all_families"].append(f"{fam_tuple[0]}{num}")

    fam_lists["plant_cell_wall"] = ["GH5", "GH6", "GH7", "GH8", "GH9", "GH10", "GH11", "GH12", "GH26", "GH28", "GH44",
                                    "GH45", "GH48", "GH53", "GH88", "GH95", "GH16", "GH17", "GH74", "GH81", "GH23",
                                    "GH27", "GH33", "GH51", "GH54", "GH62", "GH67", "GH77", "GH78", "GH84", "GH103",
                                    "GH106", "GH146",  "GH1", "GH2", "GH3", "GH13", "GH18", "GH20", "GH27", "GH29",
                                    "GH31", "GH32", "GH35", "GH38", "GH39", "GH42", "GH43", "GH52", "GH57", "GH92",
                                    "GH127", "GH130", "GH137", "GH138", "GH139", "GH141", "GH142", "GH143", "GH147"]

    fam_lists["host-specific(base)"] = ["GH2", "GH16", "GH18", "GH20", "GH29", "GH31", "GH33", "GH36", "GH84", "GH89",
                                        "GH95", "GH98", "GH101", "GH109", "GH110", "GH112", "GH123"]

    fam_lists["host-specific(extra)"] = ["GH1", "GH3", "GH5", "GH23", "GH32", "GH34", "GH35", "GH38", "GH39", "GH42",
                                         "GH43", "GH67", "GH77", "GH83", "GH85", "GH88", "GH94", "GH129", "GH130",
                                         "GH136", "GH139", "GH141", "GH151", "GH156"]

    fam_lists["mucus"] = ["GH2", "GH16", "GH16_3", "GH16_24", "GH18", "GH20", "GH29", "GH31", "GH33", "GH34", "GH36",
                          "GH84", "GH85", "GH89", "GH95", "GH98", "GH101", "GH109", "GH110", "GH112", "GH123", "GH129",
                          "GH156"]

    fam_lists["milk_glycans"] = ["GH2", "GH16", "GH18", "GH20", "GH29", "GH31", "GH33", "GH34", "GH35", "GH36", "GH84",
                                 "GH85", "GH89", "GH95", "GH98", "GH101", "GH109", "GH110", "GH112", "GH123", "GH129",
                                 "GH141", "GH156"]

    fam_lists["xyloglucan"] = ["GH5", "GH5_4", "GH9", "GH12", "GH26", "GH44", "GH45", "GH48", "GH16", "GH16_20", "GH74",
                               "GH3", "GH43", "CBM65", "CBM75"]

    fam_lists["example_category"] = ["GH62", "GT9", "PL9", "CE1", "AA2", "CBM4"]

    return fam_lists


def get_user_categories():
    #   Check for category file and create if necessary
    if not os.path.isfile(default_fam_lists_file_path):
        print("Default family category config file not found, creating it...")
        try:
            write_family_files()
        except IOError as error:
            print("ERROR:", error.args[0])
            print("ERROR: Cannot create default family_categories config file.\n"
                  "Check that you have proper filesystem permissions.")
            sys.exit(1)

    # load existing JSON which we will get list from
    try:
        with open(default_fam_lists_file_path, 'r', encoding="utf-8") as jsonfile:
            fam_cats = json.loads(jsonfile.read())
    except IOError as error:
        print("ERROR:", error.args[0])
        print("ERROR: Cannot load data from family_categories config file.\n"
              "Check that you have proper filesystem permissions.")
        sys.exit(1)

    return fam_cats


def get_category_list(category_name):
    fam_cats = get_user_categories()

    # get category from user or default categories
    if category_name in fam_cats:
        cat_list = fam_cats[category_name]
    else:
        default_fams = get_default_family_categories()
        if category_name in default_fams:
            cat_list = default_fams[category_name]
        else:
            raise UserError(f"Family category argument \"{category_name}\" is not found in family_categories.json\n"
                            f"You can add or modify custom family categories using the saccharis.add_family_category "
                            "command.")
    matcher = Matcher()
    for fam in cat_list:
        if not matcher.valid_cazy_family(fam):
            raise UserError(f"ERROR: Invalid family argument read from family category: \"{fam}\"\n"
                            f"\tPlease edit the category to either delete or edit this into a valid family: PL*, GH*, "
                            f"GT*, CE*, or AA*, where * is a number.")

    return cat_list


def load_family_list_from_file(file):
    if type(file) == str:
        family_file = open(file, 'r')
    elif type(file) == io.TextIOWrapper:
        family_file = file
    else:
        print(type(file))
        raise UserError(f"Bad variable type '{type(file)}' passed to load_family_list_from_file")

    try:
        fam_list = json.loads(family_file.read())
        if type(fam_list) == dict:
            fam_list = list(fam_list.keys())
    except JSONDecodeError as error:
        print("ERROR:", error.args[0])
        raise UserError("ERROR: Error reading user JSON file. Format is invalid. Check for extra commas at the end of "
                        "the arrays, make sure all family names are in quotes, e.g. \"GH43\", and make sure the array "
                        "is opened and closed with square brackets. An example of correct syntx is below:\n"
                        "[\n"
                        "\"GH6\",\n"
                        "\"GT9\",\n"
                        "\"PL9\",\n"
                        "\"CE1\",\n"
                        "\"AA2\",\n"
                        "\"CBM4\"\n"
                        "]\n") from error
    finally:
        family_file.close()

    matcher = Matcher()
    for fam in fam_list:
        if not matcher.valid_cazy_family(fam):
            raise UserError(f"ERROR: Invalid family argument read from file: \"{fam}\"\n"
                            f"\tPlease input a valid family: PL*, GH*, GT*, CE*, or "
                            f"AA*, where * is a number.")

    return fam_list


def save_family_iterable_json(family_iterable, out_path):
    out_dir = os.path.dirname(out_path)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # found_file = os.path.join(out_dir, re.sub(r"\.fa.*", "_families.json", os.path.basename(input_fasta)))
    with open(out_path, 'w', encoding="utf-8") as jsonfile:
        json.dump(family_iterable, jsonfile, ensure_ascii=False, indent=4)


def write_family_files(out_folder=None, data=None):
    #     todo: fix bug: default config files are not removed when program uninstalled, database files probably aren't
    #      either. BUT WAIT, is this even a bug? Don't we want config and database files to persist through software
    #      updates?
    #
    if out_folder is None:
        out_folder = folder_config
        print(f"Writing {'default ' if data is None else ''}family files to default folder location...")

    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    all_families_file_path = os.path.join(out_folder, all_families_filename)
    fam_lists_file_path = os.path.join(out_folder, fam_lists_filename)
    if data is None:
        data = get_default_family_categories()

    print(f"Writing all CAZyme family category lists to file: {fam_lists_file_path}")
    with open(fam_lists_file_path, 'w', encoding="utf-8") as file:
        json.dump(data, file, ensure_ascii=False, indent=4)

    print(f"Writing all CAZyme family list to file: {all_families_file_path}")
    with open(all_families_file_path, 'w', encoding="utf-8") as file:
        json.dump(data["all_families"], file, ensure_ascii=False, indent=4)


def show_categories():
    for category in get_user_categories():
        print("Category name:", category)
        for family in get_category_list(category):
            print(f"\t{family}")


def cli_append_user_family():
    parser = argparse.ArgumentParser(description="SACCHARIS utility to add or override a user-defined category grouping"
                                                 " of CAZyme families to the list of built in family categories used "
                                                 "with the --family_category feature.")
    parser.add_argument("category_name", type=str, help="This will be the name of the new CAzyme category you "
                                                        "are adding.")
    parser.add_argument("families", nargs='+', type=str, help="This is a space separated list of families and/or "
                                                              "subfamilies to include in your new category.")
    parser.add_argument("--file", "-f", nargs=1, type=str, help="This is an optional")
    args = parser.parse_args()
    families = args.families
    category_name = args.category_name
    m = Matcher()
    for family in families:
        if not m.valid_cazy_family(family):
            print(f"ERROR: Invalid CAZyme family: {family}")
            print("Exiting...")
            sys.exit(3)

    if not os.path.isfile(default_fam_lists_file_path):
        print("Default family category config file not found, creating it...")
        try:
            write_family_files()
        except IOError as error:
            print("ERROR:", error.args[0])
            print("ERROR: Cannot create default family_categories config file.\n"
                  "Check that you have proper filesystem permissions.")
            sys.exit(1)

    # load existing JSON which we will append or replace category of
    try:
        with open(default_fam_lists_file_path, 'r', encoding="utf-8") as jsonfile:
            categories = json.loads(jsonfile.read())
    except IOError as error:
        print("ERROR:", error.args[0])
        print("ERROR: Cannot load data from default family_category config file.")
        sys.exit(1)

    # check if user wants to overwrite a pre-existing category
    write = True
    if category_name in categories:
        write = ask_yes_no("Category already exists! Would you like to overwrite it?",
                           f"Overwriting {category_name} category...",
                           f"Leaving existing {category_name} category...")

    # write new json file data
    if write:
        categories[category_name] = families
        try:
            with open(default_fam_lists_file_path, 'w', encoding="utf-8") as jsonfile:
                json.dump(categories, jsonfile, ensure_ascii=False, indent=4)
            print(f"New category \"{category_name}\" added to category list.")
        except IOError as error:
            print("ERROR:", error.args[0])
            print("ERROR: Cannot write data to default family_category config file.")
            sys.exit(1)
    else:
        print("INFO: No data written to file.")


def cli_main():
    parser = argparse.ArgumentParser(description="Utility to generate default family files")
    parser.add_argument("--out_folder", "-f", type=str, help="Folder to output family category json files to.",
                        default=None)
    args = parser.parse_args()

    try:
        write_family_files(args.out_folder)
    except IOError as error:
        print("ERROR:", error.args[0])
        print("Family files could not be written! Check that you have permissions to write to the path.")


if __name__ == "__main__":
    cli_append_user_family()
    # cli_show_categories()

