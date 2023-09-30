import unittest
from tests_package.CazyTests import CazyTestCase
from tests_package.HelperTests import HelperTestCase
from tests_package.PruneTests import PruneTestCase
from tests_package.UserFastaRenameTest import UserRenameFastaTestCase
from tests_package.UserInputTesting import UserInputTestCase
from tests_package.DownloadTests import DownloadTestCase
from tests_package.NCBITests import NCBITestCase
from tests_package.IntegrationTests import IntegrationTestCase
from tests_package.AAModelTests import AAModelTestCase
from tests_package.TreeTests import TreeTestCase
from tests_package.RenderingTests import RenderingTestCase
from tests_package.GUITests import GUITestCase
from tests_package.ParseTests import ParseTestCase


def saccharis_test_suite():
    suite = unittest.TestSuite()
    loader = unittest.TestLoader()
    suite.addTest(loader.loadTestsFromTestCase(CazyTestCase))
    suite.addTest(loader.loadTestsFromTestCase(NCBITestCase))
    suite.addTest(loader.loadTestsFromTestCase(DownloadTestCase))
    suite.addTest(loader.loadTestsFromTestCase(PruneTestCase))
    suite.addTest(loader.loadTestsFromTestCase(UserRenameFastaTestCase))
    suite.addTest(loader.loadTestsFromTestCase(UserInputTestCase))
    suite.addTest(loader.loadTestsFromTestCase(AAModelTestCase))
    suite.addTest(loader.loadTestsFromTestCase(HelperTestCase))
    suite.addTest(loader.loadTestsFromTestCase(IntegrationTestCase))
    suite.addTest(loader.loadTestsFromTestCase(GUITestCase))
    suite.addTest(loader.loadTestsFromTestCase(TreeTestCase))
    suite.addTest(loader.loadTestsFromTestCase(RenderingTestCase))
    suite.addTest(loader.loadTestsFromTestCase(ParseTestCase))
    return suite


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(saccharis_test_suite())
