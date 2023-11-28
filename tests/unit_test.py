import unittest

class TestSum(unittest.TestCase):
    def test_the_function(self):
        lat = 16
        lon = 45
        data = [lat, lon]
        result = function_to_be_called()
        wanted_result = ""
        self.assertEqual(result, wanted_result)
        self.assertTrue(result)


if __name__ == '__main__':
    unittest.main()