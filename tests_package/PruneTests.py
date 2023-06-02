import json
import os
import shutil
import unittest
from inspect import getsourcefile

import Bio.SeqIO

from saccharis.ExtractAndPruneCAZymes import filter_prune

tests_folder = os.path.dirname(getsourcefile(lambda: 0))
testfiles_folder = os.path.join(tests_folder, "test_files")
test_out_folder = os.path.join(tests_folder, "test_files", "temp")


class PruneTestCase(unittest.TestCase):

    def setUp(self) -> None:
        if not os.path.exists(test_out_folder):
            os.mkdir(test_out_folder)

    def tearDown(self) -> None:
        print("Deleting temp files")
        shutil.rmtree(test_out_folder)

    def helper_filter_prune(self, cazy_family, input_fasta_filename, input_hmmer_filename,
                            ref_pruned_filename, ref_bounds_filename):
        pruned, bounds = filter_prune(os.path.join(testfiles_folder, input_fasta_filename),
                                      os.path.join(testfiles_folder, input_hmmer_filename),
                                      cazy_family, test_out_folder, "dbcan",
                                      prune=True, accession_dict=None)
        prunefile = os.path.join(testfiles_folder, ref_pruned_filename)
        boundsfile = os.path.join(testfiles_folder, ref_bounds_filename)
        pruned_from_file = list(Bio.SeqIO.parse(prunefile, "fasta"))
        with open(boundsfile, 'r', encoding="utf-8") as f:
            bounds_from_file = json.loads(f.read())
        for i, item in enumerate(pruned):
            self.assertEqual(item.seq, pruned_from_file[i].seq)
        for key, value in bounds.items():
            self.assertListEqual(list(value), bounds_from_file[key])

    def test_filter_prune_gh2(self):
        self.helper_filter_prune("GH2", "GH2_filter_test.fasta", "GH2_hmmer_filter_test.out", "prunedGH2.fasta",
                                 "boundsGH2.json")

    def test_filter_prune_gh16(self):
        self.helper_filter_prune("GH16", "GH16_CHARACTERIZED_cazy_filter_test.fasta", "GH16_filter_test_hmmer.out",
                                 "prunedGH16.fasta", "boundsGH16.json")


if __name__ == '__main__':
    unittest.main()
