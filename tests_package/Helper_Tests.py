import unittest
from saccharis.utils.Formatting import format_time


class HelperTestCase(unittest.TestCase):

    def test_format_hours(self):
        self.assertEqual("*\t 1 hours, 5 minutes, 17 seconds to run", format_time(3917))

    def test_format_minutes(self):
        self.assertEqual("*\t 5 minutes, 17 seconds to run", format_time(317))

    def test_format_seconds(self):
        self.assertEqual("*\t 17.0 seconds to run", format_time(17))

    def test_rounding(self):
        self.assertEqual("*\t 4 minutes, 55 seconds to run", format_time(295))
