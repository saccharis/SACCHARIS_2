import unittest
from CazyTests import CazyTestCase
from Helper_Tests import HelperTestCase
from PruneTests import PruneTestCase
from UserFastaRenameTest import UserRenameFastaTestCase
from UserInputTesting import UserInputTestCase
from DownloadTests import DownloadTestCase
from NCBITests import NCBITestCase
from IntegrationTests import IntegrationTestCase
from AAModelTests import AAModelTestCase


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
    return suite


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(saccharis_test_suite())
