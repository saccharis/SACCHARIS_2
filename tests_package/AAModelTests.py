import os
import shutil
import unittest
from saccharis import ChooseAAModel
from saccharis.Cazy_Scrape import Mode
from inspect import getsourcefile

tests_folder = os.path.dirname(getsourcefile(lambda: 0))
testfiles_folder = os.path.join(tests_folder, "test_files")
test_out_folder = os.path.join(tests_folder, "test_files", "temp")


class AAModelTestCase(unittest.TestCase):

    def setUp(self) -> None:
        self.aligned_path = os.path.join(testfiles_folder, "GH16_CHARACTERIZED.muscle_aln.phyi")
        if not os.path.exists(test_out_folder):
            os.mkdir(test_out_folder)

    def tearDown(self) -> None:
        print("Deleting temp files")
        shutil.rmtree(test_out_folder)

    def test_fasttree(self):
        best = ChooseAAModel.compute_best_model(self.aligned_path, None, "GH16", test_out_folder, 208,
                                                Mode.CHARACTERIZED, "MF", os.cpu_count(),
                                                ChooseAAModel.TreeBuilder.FASTTREE, force_update=True, user_run=None,
                                                use_modelTest=True)
        self.assertEqual(best, "gamma-wag")

    def test_raxml(self):
        best = ChooseAAModel.compute_best_model(self.aligned_path, None, "GH16", test_out_folder, 208,
                                                Mode.CHARACTERIZED, "MF", os.cpu_count(),
                                                ChooseAAModel.TreeBuilder.RAXML, force_update=True, user_run=None,
                                                use_modelTest=True)
        self.assertEqual(best, "PROTGAMMAIWAG")

    def test_fasttree_prottest(self):
        # todo: delete this test when prottest support is removed
        best = ChooseAAModel.compute_best_model(self.aligned_path, None, "GH16", test_out_folder, 208,
                                                Mode.CHARACTERIZED, "MF", os.cpu_count(),
                                                ChooseAAModel.TreeBuilder.FASTTREE, force_update=True, user_run=None,
                                                prottest_folder=r"\\wsl$\Ubuntu\usr\local\prottest-3.4.2",
                                                use_modelTest=False)
        self.assertEqual(best, "gamma-wag")

    def test_raxml_prottest(self):
        # todo: delete this test when prottest support is removed
        best = ChooseAAModel.compute_best_model(self.aligned_path, None, "GH16", test_out_folder, 208,
                                                Mode.CHARACTERIZED, "MF", os.cpu_count(),
                                                ChooseAAModel.TreeBuilder.RAXML, force_update=True, user_run=None,
                                                use_modelTest=False)
        self.assertEqual(best, "PROTGAMMAWAG")


if __name__ == '__main__':
    unittest.main()
