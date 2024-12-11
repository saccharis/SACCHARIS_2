import unittest

from tests_package.CazyTests import CazyTestCase
from tests_package.ConfigTests import ConfigTestCase
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
    """Create and return the test suite."""
    suite = unittest.TestSuite()
    loader = unittest.TestLoader()

    test_cases = [
        CazyTestCase,
        NCBITestCase,
        DownloadTestCase,
        PruneTestCase,
        UserRenameFastaTestCase,
        UserInputTestCase,
        AAModelTestCase,
        HelperTestCase,
        IntegrationTestCase,
        GUITestCase,
        TreeTestCase,
        RenderingTestCase,
        ParseTestCase,
        ConfigTestCase
    ]

    # Add all the test cases to the suite
    for test_case in test_cases:
        suite.addTest(loader.loadTestsFromTestCase(test_case))
    return suite


if __name__ == '__main__':
    unittest.TextTestRunner().run(saccharis_test_suite())
