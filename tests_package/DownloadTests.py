import os
import shutil
import sys
import unittest
from inspect import getsourcefile
from unittest import mock

from saccharis.utils.DatabaseDownload import download_database, cli_update_hmms

tests_folder = os.path.dirname(getsourcefile(lambda: 0))
test_out_folder = os.path.join(tests_folder, "test_files", "temp")


def get_hmm_suffices(filename: str):
    return [f"{filename}.h3f", f"{filename}.h3i", f"{filename}.h3m", f"{filename}.h3p"]


class DownloadTestCase(unittest.TestCase):

    def setUp(self) -> None:
        if not os.path.exists(test_out_folder):
            os.mkdir(test_out_folder)

    def tearDown(self) -> None:
        print("Deleting temp files")
        shutil.rmtree(test_out_folder)

    @unittest.skipUnless(sys.gettrace(), "expensive test, should not be run automatically for CI/CD")
    def test_db_download(self) -> None:
        files = ["CAZy.dmnd", "PUL.faa", "tcdb.dmnd", "dbCAN-PUL.tar.gz", "dbCAN.txt"]
        files += get_hmm_suffices("dbCAN.txt")
        files += get_hmm_suffices("dbCAN_sub.hmm")
        files += get_hmm_suffices("stp.hmm")
        files += get_hmm_suffices("tf-1.hmm")
        files += get_hmm_suffices("tf-2.hmm")

        download_database(test_out_folder)
        for file in files:
            self.assertTrue(os.path.isfile(os.path.join(test_out_folder, file)))
        pass

    def test_cli_update_hmms(self) -> None:
        testargs = ["saccharis.update_db", "-k", "-v", "--default_directory"]
        with mock.patch.object(sys, 'argv', testargs):
            cli_update_hmms()

    def test_cli_update_hmms_bad_dir_args(self) -> None:
        testargs = ["saccharis.update_db", "-k", "-v", "-d", "~/fake/dir", "--default_directory"]
        with mock.patch.object(sys, 'argv', testargs):
            with self.assertRaises(SystemExit):
                cli_update_hmms()
