###############################################################################
#
#
# SACCHARIS 2.0 author: Alexander Fraser, https://github.com/AlexSCFraser
# License: GPL v3
###############################################################################
# Built in libraries
import argparse
import importlib.metadata
import inspect
import json
import os
import pathlib
import shutil
import time
from importlib.metadata import version
from json import JSONDecodeError
from logging import Logger, getLogger
from subprocess import run, CalledProcessError
from dotenv import load_dotenv
import textwrap as _textwrap
# Internal imports
from saccharis.utils.UserInput import ask_yes_no


def get_version():
    try:
        return version("SACCHARIS")
    except importlib.metadata.PackageNotFoundError:
        return "dev-build"


def get_default_package_settings(logger: Logger = getLogger()):
    default_package_settings = {"user_data_folder": os.path.expanduser('~')}
    msg = f"Default package settings: {default_package_settings}"
    logger.info(msg)
    return default_package_settings


def load_from_file(path: str | os.PathLike):
    with open(path, 'r', encoding="utf-8") as sfile:
        settings = json.loads(sfile.read())
    return settings


def get_package_settings(logger: Logger = getLogger()):
    base_dir = pathlib.PurePath(inspect.getsourcefile(lambda: 0)).parents[1]
    msg = f"Base saccharis install directory: {base_dir}"
    logger.debug(msg)
    data_folder = base_dir / "data"
    config_path = data_folder / "config.json"
    if not os.path.isfile(config_path):
        msg = f"No file found at {config_path}, loading default package settings."
        logger.info(msg)
        return get_default_package_settings(logger)
    else:
        try:
            msg = f"Attempting to load package settings from {config_path}"
            logger.debug(msg)
            current_package_settings = load_from_file(config_path)
            logger.debug(f"Package settings: {current_package_settings}")
            return current_package_settings
        except JSONDecodeError:
            msg = f"JSON decode error encountered trying to load package settings from {config_path}, using defaults " \
                  f"instead."
            logger.exception(msg)
            return get_default_package_settings(logger)


def save_package_settings(new_package_settings):
    base_dir = pathlib.PurePath(inspect.getsourcefile(lambda: 0)).parents[1]
    data_folder = base_dir / "data"
    config_path = data_folder / "config.json"
    if os.path.isfile(config_path):
        try:
            old_package_settings = load_from_file(config_path)
        except JSONDecodeError:
            old_package_settings = {}
    else:
        old_package_settings = {}

    for key in old_package_settings:
        if key not in new_package_settings:
            new_package_settings[key] = old_package_settings[key]
    save_to_file(new_package_settings, config_path)


# IMPROVEMENT TODO: Consider refactoring to a global settings object with all the settings and get/set methods.
package_settings = get_package_settings()
user_dir = package_settings["user_data_folder"] if "user_data_folder" in package_settings else os.path.expanduser('~')
folder_saccharis_user = os.path.join(user_dir, "saccharis")
default_folder_config = os.path.join(folder_saccharis_user, "config")
default_settings_path = os.path.join(default_folder_config, "advanced_settings.json")
default_db_dir = os.path.join(folder_saccharis_user, "db")
folder_db = package_settings["database_folder"] if "database_folder" in package_settings else default_db_dir
folder_logs = os.path.join(folder_saccharis_user, "logs")
folder_default_output = os.path.join(folder_saccharis_user, "output")
folder_ncbi = os.path.join(folder_saccharis_user, "ncbi_downloads")


def get_db_folder(logger: Logger = getLogger()):
    msg = f"Current database_folder: {folder_db}"
    logger.debug(msg)
    return folder_db


def get_log_folder():
    return folder_logs


def get_output_folder():
    return folder_default_output


def get_config_folder():
    return default_folder_config


def get_ncbi_folder():
    return folder_ncbi


class MultilineFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.fill(paragraph, width, initial_indent=indent, subsequent_indent=indent) + \
                                  '\n'
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text


def load_from_env(gui_object=None, ask_method=None, get_method=None, show_user_method=print, skip_ask=False):
    if not os.path.isdir(default_folder_config):
        os.mkdir(default_folder_config)
    load_dotenv(os.path.join(default_folder_config, ".env"), override=True)
    if "API_KEY" in os.environ:
        api_key = os.environ['API_KEY']
    elif not skip_ask:
        show_user_method("WARNING: NCBI API KEY NOT FOUND! This will reduce speed of NCBI querying.\n"
                         "You can get an NCBI API key from your NCBI account settings page, see "
                         "https://www.ncbi.nlm.nih.gov/account/settings/")
        no_msg = "NCBI Query speed will be reduced without an API key!"
        if gui_object and ask_method:
            ask_key = ask_method("Would you like to add your NCBI API key now? (y/n):", None, no_msg, gui_object)
        else:
            ask_key = ask_yes_no("Would you like to add your NCBI API key now? (y/n):", None, no_msg)

        if ask_key:
            if gui_object and get_method:
                api_key = get_method("API key", "Please enter your NCBI API key with no spaces: ",  gui_object)
            else:
                api_key = input("Please enter your NCBI API key with no spaces: ")
            if api_key:  # check that user entered an api_key. the gui get_method returns a None type on cancel click
                api_key = api_key.strip()
                with open(os.path.join(default_folder_config, ".env"), 'a') as f:
                    f.write(f"API_KEY={api_key}\n")
                show_user_method("Note: API key has been stored in the .env file in the SACCHARIS config directory.")
                time.sleep(0.5)  # sleep the thread briefly because for some reason the email won't be written properly
                # sometimes and I think it's some kind of problem with the file not being fully closed in the OS yet.
        else:
            api_key = None
    else:
        api_key = None

    ncbi_tool = "saccharis2"

    if "EMAIL" in os.environ:
        ncbi_email = os.environ['EMAIL']
    else:
        ncbi_email = "alexscf@msl.ubc.ca"
    # note: below is old code when I thought ncbi email was per-user, but ncbi email in requests is supposed to go to
    # the software developer.
    # elif not skip_ask:
    #     show_user_method("WARNING: NCBI EMAIL NOT FOUND! This is recommended.\n"
    #                      "This email address should match the account your API key is from.")
    #     no_email_msg = "No email address added."
    #     if gui_object and ask_method:
    #         ask_email = ask_method("Would you like to add your NCBI email? (y/n):", None, no_email_msg, gui_object)
    #     else:
    #         ask_email = ask_yes_no("Would you like to add your NCBI email? (y/n):", None, no_email_msg)
    #
    #     if ask_email:
    #         if gui_object and get_method:
    #             ncbi_email = get_method("Please enter your NCBI email with no spaces.", "Email:", gui_object)
    #         else:
    #             ncbi_email = input("Please enter your NCBI email with no spaces: ")
    #         if ncbi_email:  # check that user entered an ncbi_email. the gui get_method returns None on cancel click
    #             ncbi_email = ncbi_email.strip()
    #             with open(os.path.join(folder_config, ".env"), 'a') as f:
    #                 f.write(f"EMAIL={ncbi_email}\n")
    #             show_user_method("Note: NCBI email stored in the .env file in the SACCHARIS config directory.")
    #     else:
    #         ncbi_email = None
    # else:
    #     ncbi_email = None

    return api_key, ncbi_email, ncbi_tool


def get_user_settings():
    if not os.path.isfile(default_settings_path):
        save_to_file(get_default_settings())

    try:
        user_settings = load_from_file(default_settings_path)
    except JSONDecodeError:
        user_settings = {}
    settings = get_default_settings()
    for item in settings.keys():
        if item in user_settings:
            settings[item] = user_settings[item]

    if not validate_settings(settings):
        raise UserWarning("User settings loaded from file were not valid.")

    return settings


def get_default_settings():
    settings = {"hmm_eval": 1e-15,
                "hmm_cov": 0.35,
                "genbank_query_size": 350,
                "raxml_command": "raxmlHPC-PTHREADS-AVX2",
                }

    return settings


def validate_settings(settings):

    if not type(settings["hmm_eval"]) == float:
        raise UserWarning("hmm_eval is not a float value!")
    if not type(settings["hmm_cov"]) == float:
        raise UserWarning("hmm_cov is not a float value!")
    if not type(settings["genbank_query_size"]) == int:
        raise UserWarning("genbank_query_size is not a float value!")
    if settings["genbank_query_size"] > 500:
        raise UserWarning("genbank_query_size max value is 500!")

    try:
        # if sys.platform.startswith("win"):
        #     try:
        #         rax_info = run(["wsl", settings["raxml_command"], "-v"], capture_output=True, check=True)
        #     except FileNotFoundError:
        #         raise UserWarning("Windows subsystem for linux is not installed, please install it!")
        #     except CalledProcessError:
        #         raise UserWarning(f"{settings['raxml_command']} command not found on WSL. \n"
        #                           f"raxml is not installed on your WSL installation, please install raxml on WSL with"
        #                           f" the default executable names, such as the one above.")
        # else:
        rax_info = run([settings["raxml_command"], "-v"], capture_output=True, check=True)
        if not rax_info.stdout.__contains__(b"RAxML"):
            raise UserWarning(f"Command \"{settings['raxml_command']}\" does not appear to be RAxML. Check that RAxML "
                              f"is available on path with this exact spelling. "
                              # f"Alternatively, specify the full length "
                              # f"executable path as the RAxML argument like this: "
                              # f"\n\t\"--raxml </full/path/to/RAxML_executable_file>\""
                              )

    except (FileNotFoundError, CalledProcessError) as file_e:
        raise UserWarning(f"ERROR: cannot find command \"{settings['raxml_command']}\" on path. Check that RAxML is "
                          f"available on path with this exact spelling. Alternately, if this is the wrong name for "
                          f"your raxml command, consider changing it with "
                          # f"Alternatively, specify the full length "
                          # f"executable path as the RAxML argument like this: "
                          # f"\n\t\"--raxml </full/path/to/RAxML_executable_file>\""
                          ) from file_e

    # if not (settings["raxml_command"] == "raxmlHPC-PTHREADS-AVX2" or
    #         settings["raxml_command"] == "raxmlHPC-PTHREADS-SSE3" or
    #         settings["raxml_command"] == "raxmlHPC-PTHREADS" or
    #         settings["raxml_command"] == "raxmlHPC-AVX2" or
    #         settings["raxml_command"] == "raxmlHPC-SSE3" or
    #         settings["raxml_command"] == "raxmlHPC"
    #         ):
    #     raise UserWarning("Bad raxml command in settings file")

    return True


def save_to_file(settings_dict, settings_path: str | os.PathLike = default_settings_path, api_key=None,
                 folder_config=default_folder_config):
    if not os.path.exists(os.path.dirname(settings_path)):
        os.makedirs(os.path.dirname(settings_path))

    with open(settings_path, 'w', encoding="utf-8") as sfile:
        json.dump(settings_dict, sfile, ensure_ascii=False, indent=4)

    if api_key:
        try:
            with open(os.path.join(folder_config, ".env"), 'r', encoding="utf-8") as env_file:
                lines = env_file.readlines()
        except IOError:
            lines = []

        new_lines = [(line + '\n' if line[-1] != '\n' else line) for line in lines]

        # if email:
        #     new_lines = [line for line in new_lines if line[0:5] != "EMAIL"]
        #     new_lines.append(f"EMAIL={email}\n")

        if api_key:
            new_lines = [line for line in new_lines if line[0:7] != "API_KEY"]
            new_lines.append(f"API_KEY={api_key}\n")

        with open(os.path.join(folder_config, ".env"), 'w', encoding="utf-8") as env_file:
            env_file.writelines(new_lines)


def cli_config():
    parser = argparse.ArgumentParser(description="A utility to configure advanced software_settings. If you need "
                                                 "deeper explanations of a setting, refer to the documentation of the "
                                                 "underlying tool, if any. Any software_settings changed here are "
                                                 "stored for use in future runs in the config folder: "
                                                 "~/saccharis/config",
                                     formatter_class=MultilineFormatter)
    parser.add_argument("--show", action="store_true", help="Show software_settings in console after updating values.")
    parser.add_argument('--hmm_eval', default=None, type=float, help='HMMER E Value for dbcan2')
    parser.add_argument('--hmm_cov', default=None, type=float, help='HMMER Coverage val for dbcan2')
    parser.add_argument("--querysize", "-q", default=None, type=int, help="Number of accession numbers to query genbank"
                        "with per HTTP request. Going above 350 often causes URL resolution errors. The program "
                        "automatically lowers this value when timeouts occur. It should not be necessary to specify"
                        " this value, but is left as an option in case of persistent NCBI errors.")
    parser.add_argument('--ncbi_api_key', default=None, type=str, help='This is your API key for accessing the NCBI '
                        'Entrez system. "You can get an NCBI API key from your NCBI account software_settings page, '
                        'see https://www.ncbi.nlm.nih.gov/account/settings/')
    parser.add_argument('--ncbi_email', default=None, type=str, help='This email is used for identifying your queries '
                                                                     'to NCBI in case there are errors or problems. '
                                                                     'This email address should match the account your '
                                                                     'API key is from.')
    parser.add_argument("--raxml_command", default=None, type=str, help="This is the name of the raxml command that is"
                        "available on your PATH variable. There are 6 different raxml programs that you can run, which "
                        "one you choose will be based off of the age and type of your CPU. Older "
                        "hardware may not support AVX2 or SSE3, or "
                        "multithreading. Select an option with "
                        "PTHREADS for multithreading, and AVX2 or SSE3 "
                        "for those instruction sets. AVX2 is the "
                        "fastest, SSE3 is faster than no special "
                        "instructions, and neither is slow but "
                        "supports the oldest hardware.",
                        choices=["raxmlHPC-PTHREADS-AVX2", "raxmlHPC-PTHREADS-SSE3", "raxmlHPC-PTHREADS",
                                 "raxmlHPC-AVX2", "raxmlHPC-SSE3", "raxmlHPC"])
    parser.add_argument("--user_data_folder", default=None, type=str, help="This is the folder to store SACCHARIS "
                                                                           "files containing software_settings, "
                                                                           "databases, temporary files, etc. By "
                                                                           "default, this is your user home folder "
                                                                           "as determined by expanding '~' on your "
                                                                           "operating system. "
                                                                           "For example:\n"
                                                                           "Windows: C:\\Users\\username\\ \n"
                                                                           "Linux: /home/username/\n"
                                                                           "MacOS: /Users/username/\n"
                                                                           "A folder named 'saccharis' will be created "
                                                                           "in this folder and used to store data. "
                                                                           "When changed, existing data will be "
                                                                           "moved if present."
                        )
    # old argument from Pipeline.py
    # parser.add_argument("--raxml", "-r", type=str, default="raxmlHPC-PTHREADS-AVX2", help="There are 3 different "
    #                     "raxml programs that you can run, which one you choose will be based off of the age and "
    #                     "type of your CPU. Note that this option does not do anything if your don't specify RAxML "
    #                     "via \"--tree raxml\"\n"
    #                     "\t-> Older/Slower CPU           -- raxml\n"
    #                     "\t-> Somewhat Older/Faster CPU  -- raxmlHPC-PTHREADS-SSE3\n"
    #                     "\t-> Newest Processors           -- raxmlHPC-PTHREADS-AVX2 (default)")
    parser.add_argument("--restore_defaults", "-d", action="store_true", help="If this flag is set, all values will be "
                                                                              "reset to defaults, regardless of other"
                                                                              "flags. This will NOT reset NCBI email, "
                                                                              "API key, or database folder location.")

    args = parser.parse_args()
    changes = False
    software_settings = get_user_settings()
    current_package_settings = get_package_settings()
    settings_path = default_settings_path
    config_folder = default_folder_config

    # Note: Using both expanduser and abspath so that relative paths with '.' for current working directory
    # or '~' for user's home folders are saved correctly. Otherwise, ambiguous locations might try to save/load
    # database files form incorrect locations in the future and cause database installation/update issues.
    user_data_folder = os.path.abspath(os.path.expanduser(args.user_data_folder))

    if args.hmm_eval:
        software_settings["hmm_eval"] = args.hmm_eval
        changes = True
    if args.hmm_cov:
        software_settings["hmm_cov"] = args.hmm_cov
        changes = True
    if args.querysize:
        software_settings["genbank_query_size"] = args.querysize
        changes = True
    if args.raxml_command:
        software_settings["raxml_command"] = args.raxml_command
        changes = True
    if args.user_data_folder:
        if not os.path.exists(user_data_folder):
            print(f"{user_data_folder} folder for SACCHARIS files does not exist!\nPlease choose an existing "
                  f"folder to store SACCHARIS data folder in.")
        else:
            old_user_data_folder = current_package_settings["user_data_folder"]
            old_saccharis_folder = os.path.join(old_user_data_folder, "saccharis")
            current_package_settings["user_data_folder"] = user_data_folder
            new_user_dir = current_package_settings["user_data_folder"]
            new_saccharis_folder = os.path.join(new_user_dir, "saccharis")
            config_folder = os.path.join(new_saccharis_folder, "config")
            settings_path = os.path.join(config_folder, "advanced_settings.json")

            if old_user_data_folder != current_package_settings["user_data_folder"]:
                changes = True
                shutil.copytree(old_saccharis_folder, new_saccharis_folder, dirs_exist_ok=True)
                shutil.rmtree(old_saccharis_folder)
                if current_package_settings["database_folder"] == os.path.join(old_saccharis_folder, "db"):
                    current_package_settings["database_folder"] = os.path.join(new_saccharis_folder, "db")
                save_package_settings(current_package_settings)

    if args.restore_defaults:
        software_settings = get_default_settings()
        changes = True

    try:
        validate_settings(software_settings)
        save_to_file(software_settings, settings_path=settings_path, folder_config=config_folder,
                     api_key=args.ncbi_api_key)
        if args.show:
            print(f"NCBI email: {args.ncbi_email}")
            print(f"NCBI API key: {args.ncbi_api_key}")
            print("Settings: ", json.dumps(software_settings, indent=4))
        if changes:
            print("Successfully updated advanced software_settings!")
    except UserWarning as error:
        print("ERROR", error.args[0])
        print("ERROR: Invalid software_settings, did not save changes.")


if __name__ == "__main__":
    cli_config()
