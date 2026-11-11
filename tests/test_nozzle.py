"""Check nozzle.py against textbook isentropic values (Anderson, Appendix A, gamma = 1.4)."""
import math
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "nozzle"))
import nozzle  # noqa: E402

G = 1.4


class IsentropicTables(unittest.TestCase):
    def test_area_ratio_sonic_is_one(self):
        self.assertAlmostEqual(nozzle.area_ratio(1.0, G), 1.0, places=12)

    def test_area_ratio_mach_2(self):
        self.assertAlmostEqual(nozzle.area_ratio(2.0, G), 1.6875, places=4)

    def test_pressure_ratio_mach_2(self):
        self.assertAlmostEqual(nozzle.p0_over_p(2.0, G), 7.824, places=3)

    def test_mach_from_area_inverts_area_ratio(self):
        for M in (1.2, 2.0, 3.5, 6.0):
            self.assertAlmostEqual(nozzle.mach_from_area(nozzle.area_ratio(M, G), G), M, places=6)

    def test_mach_from_pressure_inverts_p0_over_p(self):
        for M in (0.5, 1.0, 2.5):
            self.assertAlmostEqual(nozzle.mach_from_pressure(nozzle.p0_over_p(M, G), 1.0, G), M, places=9)


class Geometry(unittest.TestCase):
    def test_sharp_cone_length(self):
        # 45 deg cone: length equals radius change
        self.assertAlmostEqual(nozzle.cone_length(1.0, 0.5, 45), 0.5, places=12)

    def test_no_rounding_reduces_to_cone(self):
        self.assertAlmostEqual(nozzle.converging_length(1.0, 0.5, 30, k_round=0),
                               nozzle.cone_length(1.0, 0.5, 30), places=12)

    def test_contour_file_matches_design(self):
        path = os.path.join(os.path.dirname(__file__), "..", "nozzle_xyz.txt")
        rt, re, _ = nozzle.read_contour(path)
        self.assertAlmostEqual(2 * rt, 1.5, places=3)
        self.assertGreater(re, rt)


if __name__ == "__main__":
    unittest.main()
