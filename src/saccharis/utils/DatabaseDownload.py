import argparse
import os
import pathlib
import shutil
import subprocess
import sys

import wget
from saccharis.utils.PipelineErrors import PipelineException

from saccharis.utils.AdvancedConfig import get_db_folder, get_package_settings, save_package_settings
from saccharis.utils.Formatting import convert_path_wsl

links_last_updated = "September, 2023"
urls_and_process_and_rename = \
    [("https://bcb.unl.edu/dbCAN2/download/Databases/PUL.faa", "makeblastdb", None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL.tar.gz", "tar", None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/fam-substrate-mapping-08252022.tsv", None, None),
     ("http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_07-01-2022.xlsx", None, None),
     ("http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_07-01-2022.txt", None, None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN_sub.hmm", "hmmpress", None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/V12/tcdb.fa", "diamond", None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/V12/tf-1.hmm", "hmmpress", None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/V12/tf-2.hmm", "hmmpress", None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/V12/stp.hmm", "hmmpress", None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/V12/dbCAN-HMMdb-V12.txt", "hmmpress", "dbCAN.txt"),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/V12/CAZyDB.07262023.fa", "diamond", "CAZy.fa")
     ]

files_to_skip_deletion = ["dbCAN.txt"]


def download_and_process(url, output_folder: str | os.PathLike, process: str = None, new_filename: str = None,
                         force_download: bool = False) -> bool:
    """
    Downloads and formats a database file used for dbCAN HMMer analysis.

    :param url: url to download database file from.
    :param output_folder: folder to output downloaded processed files to
    :param process: type of process to unpack the database.
    :param new_filename: What to rename a downloaded file to prior to processing, if necessary.
    :param force_download: Forces a new download and process operation even if the files already exist.
    :return: Returns whether a file was downloaded from the URL.
    """
    downloaded = False
    if new_filename:
        output_path = os.path.join(output_folder, new_filename)
        downloaded_file = os.path.join(output_folder, os.path.basename(url))
    else:
        output_path = os.path.join(output_folder, os.path.basename(url))
        downloaded_file = output_path

    processed_filepath = f"{output_path}.h3f" if process == "hmmpress" else \
        f"{pathlib.PurePath(output_path).parent / pathlib.PurePath(output_path).stem}.dmnd" if process == "diamond" else \
            output_path

    if not os.path.exists(processed_filepath) or force_download:
        print(f"dbCAN file {processed_filepath} not found, downloading...")
        try:
            wget.download(url, output_folder)
        except TimeoutError as err:
            msg = f"Failed to download file {processed_filepath} from url {url} due to timeout."
            raise PipelineException(msg) from err
        except Exception as err:
            msg = f"Failed to download file {processed_filepath} from url {url}."
            raise PipelineException(msg) from err

        downloaded = True
        if new_filename:
            shutil.move(downloaded_file, output_path)
        assert(os.path.isfile(output_path))
        if process == "hmmpress":
            if sys.platform.startswith("win"):
                win_hmmpress_path = convert_path_wsl(output_path)
                subprocess.run(["wsl", "hmmpress", win_hmmpress_path], check=True)
            else:
                subprocess.run(["hmmpress", output_path], check=True)
            if os.path.basename(output_path) not in files_to_skip_deletion:
                os.remove(output_path)
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

    return downloaded


def cli_update_hmms():
    parser = argparse.ArgumentParser(description="Downloads the database files for HMM analysis, overwriting old files "
                                                 "if present. Note that the software may have be updated via bioconda "
                                                 "to update the internal links to download files from. Database "
                                                 f"download links last updated in {links_last_updated}")
    # fasta_filepath, bounds_file, family, output_folder, source
    parser.add_argument("--keep_existing", "-k", action="store_false", help="Use this option if you don't want newly "
                                                                            "downloaded files to replace existing "
                                                                            "files.")
    parser.add_argument("--directory", "-d", type=str, help="Use this option to download the database files into a "
                                                            "non-default directory.", default=None)

    args = parser.parse_args()
    if args.directory is None:
        old_dir = None
        dir = get_db_folder()
    else:
        old_dir = get_db_folder()
        dir = args.directory

        shutil.move(old_dir, dir)

        package_settings = get_package_settings()
        package_settings["database_folder"] = args.directory
        save_package_settings(package_settings)

    download_database(db_install_folder=dir, force_download=args.keep_existing)


def download_database(db_install_folder: str | os.PathLike = get_db_folder(), force_download: bool = False) -> int:
    """
    Downloads and formats all database files used for dbCAN HMMer analysis.

    :param db_install_folder: Folder to output downloaded and processed files to.
    :param force_download: Forces a new download and process operation even if the files already exist.
    :return: Returns the number of files downloaded and processed. The number of files created can be greater, as some
    downloaded files are processed into multiple output files.
    """

    downloaded: int = 0
    # set up folder and download dbCAN2 database files if not already present
    if not os.path.isdir(db_install_folder):
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
                                      force_download=force_download)
        if result:
            downloaded += 1

    return downloaded
