import argparse
import os
import pathlib
import time

import requests
import shutil
import subprocess
import sys
from logging import Logger, getLogger
from requests.exceptions import RequestException, ConnectionError, HTTPError
from pathlib import PurePath
from urllib3.exceptions import ReadTimeoutError

from saccharis.utils.NetworkingHelpers import resolve_hostname, get_dns_servers
from saccharis.utils.PipelineErrors import PipelineException, FileError, make_logger
from saccharis.utils.AdvancedConfig import get_db_folder, get_package_settings, save_package_settings, get_log_folder
from saccharis.utils.Formatting import convert_path_wsl

MAX_RETRIES = 10
DELAY = 30
CHUNK_SIZE: int = 4096

links_last_updated = "November, 2024"
urls_and_process_and_rename = \
    [("https://bcb.unl.edu/dbCAN2/download/Databases/V13/dbCAN-HMMdb-V13.txt", "hmmpress", "dbCAN.txt"),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/PUL.faa", "makeblastdb", None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL.tar.gz", "tar", None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/fam-substrate-mapping-08252022.tsv", None, None),
     ("http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_07-01-2022.xlsx", None, None),
     ("http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_07-01-2022.txt", None, None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN_sub.hmm", "hmmpress", None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/V12/tcdb.fa", "diamond", None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/V12/tf-1.hmm", "hmmpress", None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/V12/tf-2.hmm", "hmmpress", None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/V12/stp.hmm", "hmmpress", None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/V13/CAZyDB.07142024.fa", "diamond", "CAZy.fa")
     ]

files_to_skip_deletion = ["dbCAN.txt"]
dbcan_txt_files = ["dbCAN.txt.h3f", "dbCAN.txt.h3i", "dbCAN.txt.h3m", "dbCAN.txt.h3p"]


def download_and_process(url, output_folder: str | os.PathLike, process: str = None, new_filename: str = None,
                         force_download: bool = False, logger: Logger = getLogger()) -> bool:
    """
    Downloads and formats a database file used for dbCAN HMMer analysis.

    :param url: url to download database file from.
    :param output_folder: folder to output downloaded processed files to
    :param process: type of process to unpack the database.
    :param new_filename: What to rename a downloaded file to prior to processing, if necessary.
    :param force_download: Forces a new download and process operation even if the files already exist.
    :return: Returns whether a file was downloaded from the URL.
    """
    logger.debug(f"download_and_process() called with url:{url}; output_folder:{output_folder}; process:{process}; "
                 f"new_filename:{new_filename}; force_download:{force_download}")

    downloaded = False
    if new_filename:
        output_path = os.path.join(output_folder, new_filename)
    else:
        output_path = os.path.join(output_folder, os.path.basename(url))

    # processed_filepath = f"{output_path}.h3f" if process == "hmmpress" else \
    #     f"{PurePath(output_path).parent / PurePath(output_path).stem}.dmnd" if process == "diamond" else output_path
    processed_filepath = f"{output_path}.h3f" if process == "hmmpress" else \
        PurePath(output_path).with_suffix(".dmnd") if process == "diamond" else output_path

    if not os.path.exists(processed_filepath) or force_download:
        print(f"dbCAN file {processed_filepath} not found, downloading...")
        logger.info(f"dbCAN file {processed_filepath} not found, downloading...")

        # todo: remove this old wget file download code once confirmed requests download works
        # try:
        #     wget.download(url, output_folder)
        # except TimeoutError as err:
        #     msg = f"Failed to download file {processed_filepath} from url {url} due to timeout."
        #     logger.exception(msg)
        #     raise PipelineException(msg) from err
        # except Exception as err:
        #     msg = f"Failed to download file {processed_filepath} from url {url}."
        #     logger.debug(err)
        #     logger.exception(msg)
        #     raise PipelineException(msg) from err
        global DELAY
        global CHUNK_SIZE

        for attempt in range(MAX_RETRIES):
            try:
                response = requests.get(url, stream=True, timeout=120)
                response.raise_for_status()
                with open(output_path, 'wb') as f:
                    # consider lowering chunk size to solve proteus download issues???
                    for chunk in response.iter_content(chunk_size=CHUNK_SIZE):
                        f.write(chunk)
                logger.info(f"requests did not error on {url}")
                break  # exit retry loop on success
            except (TimeoutError, ReadTimeoutError) as err:
                if attempt < MAX_RETRIES - 1:
                    time.sleep(DELAY)
                    DELAY *= 2
                    CHUNK_SIZE = max(int(CHUNK_SIZE / 2), 512)
                    msg = f"Failed to download from {url} on attempt {attempt} due to ConnectionError. Retrying..."
                    logger.debug(msg)
                    continue
                msg = f"Failed to download file {output_path} from url {url} due to timeout. Please try again later."
                logger.debug(err)
                logger.exception(msg)
                raise PipelineException(msg) from err
            except HTTPError as err:
                msg = f"Failed to download file {output_path} from url {url} due to HTTP error. " \
                      f"Please try again later."
                logger.debug(err)
                logger.debug(resolve_hostname(url))
                logger.debug(get_dns_servers())
                logger.exception(msg)
                raise PipelineException(msg) from err
            except ConnectionError as err:
                if attempt < MAX_RETRIES - 1:
                    time.sleep(DELAY)
                    DELAY *= 2
                    CHUNK_SIZE = max(int(CHUNK_SIZE / 2), 512)
                    msg = f"Failed to download from {url} on attempt {attempt} due to ConnectionError. Retrying..."
                    logger.debug(msg)
                    continue
                msg = f"Failed to download file from url {url} due to a connection error. May be caused by DNS issues, " \
                      f"check that your system can communicate with DNS servers, sometimes VPN software or enterprise " \
                      f"network settings interfere with DNS resolution."
                logger.debug(err)
                logger.debug(resolve_hostname(url))
                logger.debug(get_dns_servers())
                logger.exception(msg)
                raise PipelineException(msg) from err
            except RequestException as err:
                msg = f"Failed to download file from url {url} due to a request error. There may be issues with the " \
                      f"server, please try again later."
                logger.debug(err)
                logger.exception(msg)
                raise PipelineException(msg) from err
            except OSError as err:
                msg = f"File I/O error occurred while accessing {output_path}. Check that you have write permissions for " \
                      f"this directory. You can change the database installation folder with the saccharis.update_db " \
                      f"command if necessary."
                logger.debug(err)
                logger.exception(msg)
                raise PipelineException(msg) from err
            except Exception as err:
                msg = f"Failed to download file {output_path} from url {url} due to generic Exception."
                logger.debug(err)
                logger.exception(msg)
                raise PipelineException(msg) from err

        downloaded = True

        # todo: delete this rename, unneeded with reworked requests download
        # if new_filename:
        #     shutil.move(downloaded_file, output_path)

        if not os.path.isfile(output_path):
            msg = f"{output_path} file does not exist! Error downloading and/or processing this file."
            logger.error(msg)
            raise FileError(msg)

        if process == "hmmpress":
            if sys.platform.startswith("win"):
                win_hmmpress_path = convert_path_wsl(output_path)
                subprocess.run(["wsl", "hmmpress", "-f", win_hmmpress_path], check=True)
            else:
                subprocess.run(["hmmpress", "-f", output_path], check=True)
            if os.path.basename(output_path) not in files_to_skip_deletion:
                os.remove(output_path)
                logger.debug(f"Removed {output_path}")
        elif process == "tar":
            if sys.platform.startswith("win"):
                win_tar_path = convert_path_wsl(output_path)
                subprocess.run(["wsl", "tar", "xvf", win_tar_path], check=True)
            else:
                subprocess.run(["tar", "xvf", output_path], check=True)
        elif process == "makeblastdb":
            subprocess.run(["makeblastdb", "-in", output_path, "-dbtype", "prot"], check=True)
        elif process == "diamond":
            output_path_obj = pathlib.PurePath(output_path)
            diamond_output_path = output_path_obj.parent / output_path_obj.stem
            subprocess.run(["diamond", "makedb", "--in", output_path, "-d", diamond_output_path], check=True)
            if os.path.basename(output_path) not in files_to_skip_deletion:
                os.remove(output_path)
                logger.debug(f"Removed {output_path}")
    else:
        logger.debug(f"Size of {processed_filepath}: {os.path.getsize(processed_filepath)}")
        logger.info(f"{processed_filepath} already exists!")

    # Below code exists to create a dbCAN.txt dummy file if it's accidentally deleted. It's not strictly
    # necessary, but run-dbcan checks for it and fails if not present even if the dbCAN.txt.h3* files are present.
    # In theory this is unneeded since we don't delete files in the files_to_skip_deletion list, but some systems
    # delete dbCAN.txt anyway for mysterious reasons, so I added this check.
    if new_filename == "dbCAN.txt" and not os.path.isfile(output_path):
        dbcan_txt_exists = [os.path.exists(os.path.join(output_folder, dbcan_txt_file))
                                for dbcan_txt_file in dbcan_txt_files]
        if all(dbcan_txt_exists):
            pathlib.Path(output_path).touch(exist_ok=True)

    return downloaded


def cli_update_hmms():
    """
    Command line interface to download/update database files used for dbCAN HMMer analysis. Can also be used to
    move/install the database in a custom location, useful for systems where the home folder is limited in size
    (e.g. computing clusters). User calls this function using the 'saccharis.update_db' command.
    """
    parser = argparse.ArgumentParser(description="Downloads the database files for HMM analysis, overwriting old files "
                                                 "if present. Note that the software may have be updated via bioconda "
                                                 "to update the internal links to download files from. Database "
                                                 f"download links last updated in {links_last_updated}")
    # fasta_filepath, bounds_file, family, output_folder, source
    parser.add_argument("--keep_existing", "-k", action="store_true", help="Use this option if you don't want newly "
                                                                            "downloaded files to replace existing "
                                                                            "files.")
    directory_group = parser.add_mutually_exclusive_group(required=False)
    directory_group.add_argument("--directory", "-d", type=str, help="Use this option to download the database files "
                                                                     "into a non-default directory.", default=None)
    directory_group.add_argument("--default_directory", action="store_true", help="Use this option to download the "
                                                                                  "database files into the default "
                                                                                  "directory.")
    parser.add_argument("--verbose", "-v", action="store_true", help="Use this option to print more detailed "
                                                                     "information to the console.")
    args = parser.parse_args()

    logger = make_logger("CLILogger", get_log_folder(), "cli_logs.txt", verbose=args.verbose)

    msg = f"cli_update_hmms() called with user arguments directory:{args.directory} and " \
          f"keep_existing:{args.keep_existing}"
    logger.info(msg)

    if args.directory is None:
        old_dir = None
        dir = get_db_folder(logger)
        msg = f"User did not specify directory, will attempt to download into current saccharis db directory: {dir}"
        logger.info(msg)
    else:
        old_dir = get_db_folder(logger)
        # Note: Using both expanduser and abspath so that relative paths with '.' for current working directory
        # or '~' for user's home folders are saved correctly. Otherwise, ambiguous locations might try to save/load
        # database files form incorrect locations in the future and cause database installation/update issues.
        dir = os.path.abspath(os.path.expanduser(args.directory))

        msg = f"User specified a directory: {dir} and current saccharis db directory exists at: {dir} ." \
              f"Attempting to move existing files to new directory."
        logger.info(msg)

        shutil.move(old_dir, dir)

        package_settings = get_package_settings(logger)
        package_settings["database_folder"] = dir
        msg = f"Setting database_folder in package_settings to {dir} and attempting to save."
        logger.info(msg)
        save_package_settings(package_settings)

    try:
        logger.info(f"Attempting to download database files into {dir}")
        download_database(db_install_folder=dir, force_download=not args.keep_existing, logger=logger)
    except PipelineException:
        logger.error(PipelineException.args[0])
        logger.debug(PipelineException.__traceback__)


def download_database(db_install_folder: str | os.PathLike = get_db_folder(), force_download: bool = False,
                      logger: Logger = getLogger()) -> int:
    """
    Downloads and formats all database files used for dbCAN HMMer analysis.

    :param db_install_folder: Folder to output downloaded and processed files to.
    :param force_download: Forces a new download and process operation even if the files already exist.
    :return: Returns the number of files downloaded and processed. The number of files created can be greater, as some
    downloaded files are processed into multiple output files.
    """

    msg = f"download_database() called with db_install_folder:{db_install_folder} and force_download:{force_download}"
    logger.debug(msg)
    downloaded: int = 0
    # set up folder and download dbCAN2 database files if not already present
    if not os.path.isdir(db_install_folder):
        msg = f"{db_install_folder} not found, creating directory..."
        logger.info(msg)
        os.makedirs(db_install_folder, 0o755)

    # todo: delete below after confirmed it's been replaced properly
    # if not os.path.exists(os.path.join(db_install_folder, "CAZy.dmnd")) or force_download:
    #     print("dbCAN2 file 1/12 not found, downloading...")
    #     wget.download("http://bcb.unl.edu/dbCAN2/download/Databases/V11/CAZyDB.08062022.fa", db_install_folder)
    #     if sys.platform.startswith("win"):
    #         # todo: change this to occur when command is missing?? kind of unnecessary, since diamond is availabe on windows via manual install
    #         diamond_db_inpath = convert_path_wsl(os.path.join(db_install_folder, "CAZyDB.08062022.fa"))
    #         diamond_db_outpath = convert_path_wsl(os.path.join(db_install_folder, "CAZy"))
    #         subprocess.run(["wsl", "diamond", "makedb", "--in", diamond_db_inpath, "-d", diamond_db_outpath])
    #     else:
    #         subprocess.run(["diamond", "makedb", "--in", os.path.join(db_install_folder, "CAZyDB.08062022.fa"), "-d",
    #                         os.path.join(db_install_folder, "CAZy")])
    #     os.remove(os.path.join(db_install_folder, "CAZyDB.08062022.fa"))  # 1 gigabyte file not needed, no point keeping

    # todo: delete below after confirmed it's been replaced properly
    # if not os.path.exists(os.path.join(db_install_folder, "dbCAN.txt")) or force_download:
    #     print("dbCAN2 file 2 not found, downloading...")
    #     wget.download("https://bcb.unl.edu/dbCAN2/download/Databases/V11/dbCAN-HMMdb-V11.txt", db_install_folder)
    #     shutil.move(os.path.join(db_install_folder, "dbCAN-HMMdb-V11.txt"), os.path.join(db_install_folder, "dbCAN.txt"))
    #     if sys.platform.startswith("win"):
    #         win_hmmpress_path = subprocess.run(
    #             ["wsl", "wslpath", "'" + os.path.join(db_install_folder, "dbCAN.txt") + "'"],
    #             capture_output=True, check=True).stdout.decode().strip()
    #         subprocess.run(["wsl", "hmmpress", win_hmmpress_path])
    #     else:
    #         subprocess.run(["hmmpress", os.path.join(db_install_folder, "dbCAN.txt")])

    for url, process_type, new_file in urls_and_process_and_rename:
        result = download_and_process(url, db_install_folder, process_type, new_filename=new_file,
                                      force_download=force_download, logger=logger)
        if result:
            downloaded += 1
        logger.debug(f"{downloaded} files downloaded.")

    return downloaded
