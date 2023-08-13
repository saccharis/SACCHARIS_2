import shutil
import unittest
import os
from inspect import getsourcefile

from saccharis.Pipeline import single_pipeline
from saccharis.Cazy_Scrape import Mode

tests_folder = os.path.dirname(getsourcefile(lambda: 0))
test_out_folder = os.path.join(tests_folder, "test_files", "temp")


class IntegrationTestCase(unittest.TestCase):

    def setUp(self) -> None:
        if not os.path.exists(test_out_folder):
            os.mkdir(test_out_folder)

    def tearDown(self) -> None:
        print("Deleting temp files")
        shutil.rmtree(test_out_folder)

    def run_pipeline(self, family, scrape_mode: Mode):
        single_pipeline(family=family, output_folder=test_out_folder, scrape_mode=scrape_mode, skip_user_ask=True,
                        force_update=True, verbose=True)
        self.assertTrue(os.path.exists(os.path.join(test_out_folder, f"{family}_{scrape_mode.name}_ALL_DOMAINS", f"{family}_{scrape_mode.name}_ALL_DOMAINS.json")))
        self.assertTrue(os.path.exists(os.path.join(test_out_folder, f"{family}_{scrape_mode.name}_ALL_DOMAINS", f"{family}_{scrape_mode.name}_ALL_DOMAINS_FASTTREE.tree")))

    # def test_pipeline(self):
    #     single_pipeline(family="PL9", output_folder=test_out_folder, scrape_mode=Mode.CHARACTERIZED, skip_user_ask=True,
    #                     force_update=True, verbose=True)
    #     self.assertTrue(os.path.exists(os.path.join(test_out_folder, "PL9_CHARACTERIZED_ALL_DOMAINS", "PL9_CHARACTERIZED_ALL_DOMAINS.json")))
    #     self.assertTrue(os.path.exists(os.path.join(test_out_folder, "PL9_CHARACTERIZED_ALL_DOMAINS", "PL9_CHARACTERIZED_ALL_DOMAINS_FASTTREE.tree")))

    def test_PL9(self):
        self.run_pipeline("PL9", Mode.CHARACTERIZED)

    def test_GH5_25(self):
        self.run_pipeline("GH5_25", Mode.CHARACTERIZED)


if __name__ == '__main__':
    unittest.main()
