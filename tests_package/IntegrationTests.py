import json
import shutil
import unittest
import os
from inspect import getsourcefile

from saccharis.Pipeline import single_pipeline
from saccharis.Cazy_Scrape import Mode
from saccharis.ChooseAAModel import TreeBuilder
from saccharis.utils.Formatting import CazymeMetadataRecord
from saccharis.utils.PipelineErrors import AAModelError

tests_folder = os.path.dirname(getsourcefile(lambda: 0))
test_out_folder = os.path.join(tests_folder, "test_files", "temp")
small_user_testfile = os.path.join(tests_folder, "test_files", "user_test_GH102_UserFormat.fasta")
partial_modeltest_folder = os.path.join(tests_folder, "test_files", "partial_run_modeltest", "PL9_CHARACTERIZED_ALL_DOMAINS")
sheep3_user_testfile = os.path.join(tests_folder, "test_files", "Sheep_3_protein.fasta")
sheep4_user_testfile = os.path.join(tests_folder, "test_files", "Sheep_4_protein.fasta")


class IntegrationTestCase(unittest.TestCase):

    def setUp(self) -> None:
        if not os.path.exists(test_out_folder):
            os.mkdir(test_out_folder)

    def tearDown(self) -> None:
        print("Deleting temp files")
        shutil.rmtree(test_out_folder)

    def run_pipeline(self, family, scrape_mode: Mode, tree_program: TreeBuilder = TreeBuilder.FASTTREE,
                     user_file: str = None, force_update=True, render_trees=False):
        single_pipeline(family=family, output_folder=test_out_folder, scrape_mode=scrape_mode, skip_user_ask=True,
                        force_update=force_update, verbose=True, tree_program=tree_program, user_file=user_file,
                        render_trees=render_trees)
        domain_folder = f"{family}_{scrape_mode.name}_ALL_DOMAINS"
        file_prefix = f"{family}_{scrape_mode.name}_ALL_DOMAINS{'_UserRun00000' if user_file else ''}"
        tree_prog = tree_program.name
        json_path = os.path.join(test_out_folder, domain_folder, f"{file_prefix}.json")
        self.assertTrue(os.path.exists(json_path))
        self.assertTrue(os.path.exists(os.path.join(test_out_folder, domain_folder, f"{file_prefix}_{tree_prog}.tree")))
        with open(json_path, 'r', encoding="utf-8") as meta_json:
            cazyme_dict = json.loads(meta_json.read())
            final_metadata_dict = {id: CazymeMetadataRecord(**record) for id, record in cazyme_dict.items()}
        # asserts that there are no exactly overlapping modules from multiple genes
        for record in final_metadata_dict:
            if record.__contains__("<1>"):
                record_2 = final_metadata_dict[record].protein_id + "<2>"
                self.assertFalse(final_metadata_dict[record].module_start == final_metadata_dict[record_2].module_start)
                self.assertFalse(final_metadata_dict[record].module_end == final_metadata_dict[record_2].module_end)
        if render_trees:
            tree_files = ["Basic_circular_tree.pdf", "basic_circular_tree_bootstrap.pdf",
                          "basic_circular_with_domain.pdf", "basic_circular_domain_bootstrap.pdf",
                          "basic_circular_domain_ECNo.pdf", "basic_circular_domain_ECno_numeric.pdf",
                          "Rplots.pdf"]
            for filename in tree_files:
                self.assertTrue(os.path.isfile(os.path.join(test_out_folder, domain_folder, filename)))

    def test_PL9(self):
        self.run_pipeline("PL9", Mode.CHARACTERIZED)

    def test_PL9_raxml(self):
        self.run_pipeline("PL9", Mode.CHARACTERIZED, tree_program=TreeBuilder.RAXML, user_file=small_user_testfile)

    def test_PL9_raxml_ng(self):
        self.run_pipeline("PL9", Mode.CHARACTERIZED, tree_program=TreeBuilder.RAXML_NG, user_file=small_user_testfile,
                          render_trees=True)

    def test_PL9_raxml_structure(self):
        self.run_pipeline("PL9", Mode.STRUCTURE, tree_program=TreeBuilder.RAXML, user_file=small_user_testfile)

    def test_GH5_25(self):
        self.run_pipeline("GH5_25", Mode.CHARACTERIZED)

    def test_GH5_4(self):
        self.run_pipeline("GH5_4", Mode.CHARACTERIZED)

    def test_bad_partial_modeltest_pl9(self):
        out_folder = os.path.join(test_out_folder, "PL9_CHARACTERIZED_ALL_DOMAINS")
        if os.path.exists(out_folder):
            shutil.rmtree(out_folder)
        shutil.copytree(partial_modeltest_folder, out_folder)
        self.assertRaises(AAModelError, self.run_pipeline, "PL9", Mode.CHARACTERIZED, force_update=False)

    def test_duplicate_user_files(self):
        self.run_pipeline("GH1", Mode.CHARACTERIZED, render_trees=True, user_file=sheep3_user_testfile)
        #todo: add sheep4 user test file in dev22 version


if __name__ == '__main__':
    unittest.main()
