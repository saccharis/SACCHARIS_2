import os
import pathlib
import shutil
import subprocess
import sys

import wget

from saccharis.utils.AdvancedConfig import get_db_folder
from saccharis.utils.Formatting import convert_path_wsl

urls_and_process_and_rename = \
    [("https://bcb.unl.edu/dbCAN2/download/Databases/PUL.faa", "makeblastdb", None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL.tar.gz", "tar", None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/V11/tcdb.fa", "diamond", None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/fam-substrate-mapping-08252022.tsv", None, None),
     ("http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_07-01-2022.xlsx", None, None),
     ("http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_07-01-2022.txt", None, None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/V11/tf-1.hmm", "hmmpress", None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/V11/tf-2.hmm", "hmmpress", None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/V11/stp.hmm", "hmmpress", None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN_sub.hmm", "hmmpress", None),
     ("https://bcb.unl.edu/dbCAN2/download/Databases/V11/dbCAN-HMMdb-V11.txt", "hmmpress", "dbCAN.txt"),
     ("http://bcb.unl.edu/dbCAN2/download/Databases/V11/CAZyDB.08062022.fa", "diamond", "CAZy.fa")
     ]


def download_and_process(url, output_folder: str | os.PathLike, process: str = None, new_filename: str = None,
                         force_download: bool = False):
    if new_filename:
        output_path = os.path.join(output_folder, new_filename)
        downloaded_file = os.path.join(output_folder, os.path.basename(url))
    else:
        output_path = os.path.join(output_folder, os.path.basename(url))
        downloaded_file = output_path

    processed_filepath = f"{output_path}.h3f" if process == "hmmpress" else \
                         f"{pathlib.PurePath(output_path).stem}.dmnd" if process == "diamond" else \
                         output_path

    if not os.path.exists(processed_filepath) or force_download:
        print("dbCAN file not found, downloading...")
        wget.download(url, output_folder)
        if new_filename:
            shutil.move(downloaded_file, output_path)
        assert(os.path.isfile(output_path))
        if process == "hmmpress":
            if sys.platform.startswith("win"):
                win_hmmpress_path = convert_path_wsl(output_path)
                subprocess.run(["wsl", "hmmpress", win_hmmpress_path], check=True)
            else:
                subprocess.run(["hmmpress", output_path], check=True)
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
            os.remove(output_path)


def update_hmms():
    download_database(force_download=True)


def download_database(db_install_folder: str | os.PathLike = get_db_folder(), force_download=False):

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
        download_and_process(url, db_install_folder, process_type, new_filename=new_file, force_download=force_download)
