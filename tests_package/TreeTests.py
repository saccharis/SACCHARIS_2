import sys
import unittest
from subprocess import run


class TreeTestCase(unittest.TestCase):
    def test_raxml_ng_available(self):
        args = ["raxml-ng", "--version"]
        if sys.platform.startswith("win"):
            args.insert(0, "wsl")
        result = run(args, check=True)
        self.assertTrue(result.returncode == 0)


if __name__ == '__main__':
    unittest.main()
