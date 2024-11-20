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
    """Create and return the test suite."""
    suite = unittest.TestSuite()
    loader = unittest.TestLoader()

    # Add all the test cases to the suite
    suite.addTests(loader.loadTestsFromTestCase(CazyTestCase))
    suite.addTests(loader.loadTestsFromTestCase(NCBITestCase))
    suite.addTests(loader.loadTestsFromTestCase(DownloadTestCase))
    suite.addTests(loader.loadTestsFromTestCase(PruneTestCase))
    suite.addTests(loader.loadTestsFromTestCase(UserRenameFastaTestCase))
    suite.addTests(loader.loadTestsFromTestCase(UserInputTestCase))
    suite.addTests(loader.loadTestsFromTestCase(AAModelTestCase))
    suite.addTests(loader.loadTestsFromTestCase(HelperTestCase))
    suite.addTests(loader.loadTestsFromTestCase(IntegrationTestCase))

    return suite


if __name__ == '__main__':
    unittest.TextTestRunner().run(saccharis_test_suite())
