import os
import shutil
import sys
import unittest
from inspect import getsourcefile
from unittest import mock


from saccharis.utils.AdvancedConfig import cli_config

tests_folder = os.path.dirname(getsourcefile(lambda: 0))
test_out_folder = os.path.join(tests_folder, "test_files", "temp")


class ConfigTestCase(unittest.TestCase):

    def setUp(self) -> None:
        if not os.path.exists(test_out_folder):
            os.mkdir(test_out_folder)

    def tearDown(self) -> None:
        print("Deleting temp files")
        shutil.rmtree(test_out_folder)

    def test_config_command(self):
        args = ["saccharis.config"]
        with mock.patch.object(sys, 'argv', args):
            with self.assertRaises(SystemExit):
                cli_config()

    def test_show(self):
        args = ["saccharis.config", "--show"]
        with mock.patch.object(sys, 'argv', args):
            with self.assertRaises(SystemExit):
                cli_config()

    def test_help(self):
        args = ["saccharis.config", "--help"]
        with mock.patch.object(sys, 'argv', args):
            with self.assertRaises(SystemExit):
                cli_config()
