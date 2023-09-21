import unittest


class GUITestCase(unittest.TestCase):
    def test_UIDesign_import(self):
        from saccharis.gui import UIDesign
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
