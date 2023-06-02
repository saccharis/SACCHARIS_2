import os
import unittest
from inspect import getsourcefile

from saccharis.utils.UserFastaRename import rename_fasta_file

tests_folder = os.path.dirname(getsourcefile(lambda: 0))
testfiles_folder = os.path.join(tests_folder, "test_files")


class UserRenameFastaTestCase(unittest.TestCase):
    def test_wrong_filetype(self):
        with self.assertRaises(UserWarning) as cm:
            rename_fasta_file(os.path.join(testfiles_folder, "accession_list.csv"))
        self.assertEqual(cm.exception.args[0],
                         "File contains no valid sequences! Check that the file is in a valid FASTA format.")


if __name__ == '__main__':
    unittest.main()
