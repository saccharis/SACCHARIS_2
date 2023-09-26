import os
import shutil
import unittest
from inspect import getsourcefile

from saccharis.Rendering import render_phylogeny

tests_folder = os.path.dirname(getsourcefile(lambda: 0))
test_files_folder = os.path.join(tests_folder, "test_files")
test_out_folder = os.path.join(test_files_folder, "temp")
json_testfile = os.path.join(tests_folder, "test_files", "PL9_CHARACTERIZED_ALL_DOMAINS.json")
tree_testfile = os.path.join(tests_folder, "test_files", "PL9_CHARACTERIZED_ALL_DOMAINS_FASTTREE.tree")


class RenderingTestCase(unittest.TestCase):

    def setUp(self) -> None:
        if not os.path.exists(test_out_folder):
            os.mkdir(test_out_folder)

    def tearDown(self) -> None:
        print("Deleting temp files")
        shutil.rmtree(test_out_folder)

    def test_render(self):
        render_phylogeny(json_testfile, tree_testfile, test_out_folder)
        tree_files = ["Basic_circular_tree.pdf", "basic_circular_tree_bootstrap.pdf",
                      "basic_circular_with_domain.pdf", "basic_circular_domain_bootstrap.pdf",
                      "basic_circular_domain_ECNo.pdf", "basic_circular_domain_ECno_numeric.pdf",
                      "Rplots.pdf"]
        for filename in tree_files:
            self.assertTrue(os.path.isfile(os.path.join(test_out_folder, filename)))


if __name__ == '__main__':
    unittest.main()
