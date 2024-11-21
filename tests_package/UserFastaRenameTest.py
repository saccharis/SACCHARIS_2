import logging
import os
import shutil
import unittest
from inspect import getsourcefile
from pathlib import Path

from saccharis.ParseUserSequences import concatenate_multiple_fasta
from saccharis.utils.Formatting import rename_metadata_dict_ids
from saccharis.utils.UserFastaRename import rename_fasta_file

tests_folder = Path(os.path.dirname(getsourcefile(lambda: 0)))
testfiles_folder = tests_folder / "test_files"
test_out_folder = testfiles_folder / "temp"
sheep3_user_testfile = testfiles_folder / "Sheep_3_protein.fasta"
sheep4_user_testfile = testfiles_folder / "Sheep_4_protein.fasta"


class UserRenameFastaTestCase(unittest.TestCase):

    def setUp(self) -> None:
        logging.basicConfig(level=logging.DEBUG)
        if not os.path.exists(test_out_folder):
            os.mkdir(test_out_folder)

    def tearDown(self) -> None:
        print("Deleting temp files")
        shutil.rmtree(test_out_folder)

    def test_wrong_filetype(self):
        with self.assertRaises(UserWarning) as cm:
            rename_fasta_file(os.path.join(testfiles_folder, "accession_list.csv"))
        self.assertEqual(cm.exception.args[0],
                         "File contains no valid sequences! Check that the file is in a valid FASTA format.")

    def test_merge_files(self):
        user_path, user_merged_dict, user_seqs = concatenate_multiple_fasta([sheep3_user_testfile,
                                                                             sheep4_user_testfile],
                                                                            output_folder=test_out_folder)
        new_user_path = rename_fasta_file(user_path)
        user_merged_dict = rename_metadata_dict_ids(new_user_path, user_merged_dict)

        self.assertEqual(len(user_merged_dict), 10683)


if __name__ == '__main__':
    unittest.main()
