###############################################################################
#
#
# SACCHARIS 2.0 author: Alexander Fraser, https://github.com/AlexSCFraser
# License: GPL v3
###############################################################################
import os
import shutil
import sys
import unittest
from inspect import getsourcefile
from unittest import mock

from saccharis.ScreenUserFile import get_user_selection
from saccharis.ScreenUserFile import extract_families_hmmer
from saccharis.utils.UserInput import ask_yes_no

tests_folder = os.path.dirname(getsourcefile(lambda: 0))
testfiles_folder = os.path.join(tests_folder, "test_files")
test_out_folder = os.path.join(tests_folder, "test_files", "temp")


class UserInputTestCase(unittest.TestCase):
    def setUp(self) -> None:
        if not os.path.exists(test_out_folder):
            os.mkdir(test_out_folder)

    def tearDown(self) -> None:
        print("Deleting temp files")
        shutil.rmtree(test_out_folder)

    @mock.patch('saccharis.ScreenUserFile.input', create=True)
    def test_mock_user(self, mocked_input):
        mocked_input.side_effect = ["GH1 PL9 AA6"]
        initial_dict = {"GH1": 1, "PL9": 2, "AA6": 5}
        fam_list = get_user_selection(initial_dict)
        self.assertEqual(len(initial_dict), len(fam_list))

    @mock.patch('saccharis.ScreenUserFile.input', create=True)
    def test_category_xyloglucan(self, mocked_input):
        # initial_dict = {"GH1": 1, "PL9": 2, "AA6": 5}
        mocked_input.side_effect = ["xyloglucan"]
        initial_dict = {"GH5": 1, "PL9": 2, "GH9": 5, "GT6": 1, "PL4": 2, "AA132_21": 5, "GH19": 1, "PL78": 2}
        fam_list = get_user_selection(initial_dict)
        self.assertEqual(2, len(fam_list))

    @unittest.skipUnless(sys.gettrace(), "expensive test, should not be run automatically for CI/CD")
    def test_extract(self):
        expected_dict = {"PL9": 12, "PL9_1": 5, "PL9_4": 1}
        file_path = os.path.join(testfiles_folder, "PL9.fasta")
        # todo: make run_dbcan work on windows, then finish this test
        fam_dict = extract_families_hmmer(file_path, test_out_folder, threads=2)
        self.assertEqual(fam_dict, expected_dict)

    @mock.patch('saccharis.utils.UserInput.input', create=True)
    def test_yes_no(self, mocked_input):
        mocked_input.side_effect = ["y", "n", "invalid", "yes", "n/a", "no"]
        self.assertTrue(ask_yes_no("Input yes for test purposes:", "yes entered", "no entered"))
        self.assertFalse(ask_yes_no("Input no for test purposes:", "yes entered", "no entered"))
        self.assertTrue(ask_yes_no("Input yes for test purposes:", "yes entered", "no entered"))
        self.assertFalse(ask_yes_no("Input no for test purposes:", "yes entered", "no entered"))


if __name__ == '__main__':
    unittest.main()
