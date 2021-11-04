import unittest

from pygenutils import NumericRange, NumericRangeSet


class TestNumericRange(unittest.TestCase):

    def test_equal(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(1, 10)

        self.assertEqual(nr1, nr2)

    def test_equal_not(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(1, 11)

        self.assertNotEqual(nr1, nr2)

    def test_le(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(1, 10)

        self.assertLessEqual(nr1, nr2)
        self.assertLessEqual(nr2, nr1)

    def test_lt_end(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(1, 11)

        self.assertLess(nr1, nr2)

    def test_lt_start(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(2, 10)

        self.assertLess(nr1, nr2)

    def test_intersection(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(5, 15)
        nr = nr1 & nr2

        self.assertEqual(nr, {NumericRange(5, 10)})

    def test_intersection_edge(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(10, 15)
        nr = nr1 & nr2

        self.assertEqual(nr, {NumericRange(10, 10)})

    def test_union(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(5, 15)
        nr = nr1 | nr2

        self.assertEqual(nr, {NumericRange(1, 15)})

    def test_union_interior(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(5, 9)
        nr = nr1 | nr2

        self.assertEqual(nr, {NumericRange(1, 10)})

    def test_union_adjacent(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(11, 15)
        nr = nr1 | nr2

        self.assertEqual(nr, {NumericRange(1, 15)})

    def test_union_disjoint(self) -> None:
        nr1 = NumericRange(1, 5)
        nr2 = NumericRange(11, 15)
        nr = nr1 | nr2

        self.assertEqual(
            nr,
            {
                NumericRange(1, 5),
                NumericRange(11, 15)
            }
        )

    def test_isdisjoint(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(11, 15)

        self.assertTrue(nr1.isdisjoint(nr2))
        self.assertTrue(nr2.isdisjoint(nr1))

    def test_isdisjoint_single_overlap(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(10, 15)

        self.assertFalse(nr1.isdisjoint(nr2))
        self.assertFalse(nr2.isdisjoint(nr1))

    def test_isdisjoint_total_overlap(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(5, 9)

        self.assertFalse(nr1.isdisjoint(nr2))
        self.assertFalse(nr2.isdisjoint(nr1))

    def test_isdisjoint_separated(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(12, 15)

        self.assertTrue(nr1.isdisjoint(nr2))
        self.assertTrue(nr2.isdisjoint(nr1))

    def test_isadjacent(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(11, 15)

        self.assertTrue(nr1.isadjacent(nr2))
        self.assertTrue(nr2.isadjacent(nr1))

    def test_isadjacent_separated(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(12, 15)

        self.assertFalse(nr1.isadjacent(nr2))
        self.assertFalse(nr2.isadjacent(nr1))

    def test_isadjacent_overlap(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(9, 15)

        self.assertFalse(nr1.isadjacent(nr2))
        self.assertFalse(nr2.isadjacent(nr1))

    def test_issubset_overlap(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(9, 15)

        self.assertFalse(nr1.issubset(nr2))
        self.assertFalse(nr2.issubset(nr1))

    def test_issubset_strict(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(4, 9)

        self.assertFalse(nr1.issubset(nr2))
        self.assertTrue(nr2.issubset(nr1))

    def test_issubset_equals(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(1, 10)

        self.assertTrue(nr1.issubset(nr2))
        self.assertTrue(nr2.issubset(nr1))

    def test_issuperset_overlap(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(9, 15)

        self.assertFalse(nr1.issuperset(nr2))
        self.assertFalse(nr2.issuperset(nr1))

    def test_issuperset_strict(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(4, 9)

        self.assertTrue(nr1.issuperset(nr2))
        self.assertFalse(nr2.issuperset(nr1))

    def test_issuperset_equals(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(1, 10)

        self.assertTrue(nr1.issuperset(nr2))
        self.assertTrue(nr2.issuperset(nr1))

    def test_sub(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(9, 15)
        nr = nr1 - nr2

        self.assertEqual(nr, {NumericRange(1, 8)})

    def test_sub_disjoints(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(13, 15)
        nr = nr1 - nr2

        self.assertEqual(nr, {nr1})

    def test_sub_adjacent(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(11, 15)
        nr = nr1 - nr2

        self.assertEqual(nr, {nr1})

    def test_sub_equals(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(1, 10)
        nr = nr1 - nr2

        self.assertEqual(nr, set())

    def test_sub_left(self) -> None:
        nr1 = NumericRange(5, 10)
        nr2 = NumericRange(1, 6)
        nr = nr1 - nr2

        self.assertEqual(nr, {NumericRange(7, 10)})

    def test_sub_left_edge(self) -> None:
        nr1 = NumericRange(5, 10)
        nr2 = NumericRange(1, 10)
        nr = nr1 - nr2

        self.assertEqual(nr, set())

    def test_sub_right(self) -> None:
        nr1 = NumericRange(5, 10)
        nr2 = NumericRange(8, 15)
        nr = nr1 - nr2

        self.assertEqual(nr, {NumericRange(5, 7)})

    def test_sub_right_edge(self) -> None:
        nr1 = NumericRange(5, 10)
        nr2 = NumericRange(5, 15)
        nr = nr1 - nr2

        self.assertEqual(nr, set())

    def test_sub_superset(self) -> None:
        nr1 = NumericRange(5, 10)
        nr2 = NumericRange(1, 15)
        nr = nr1 - nr2

        self.assertEqual(nr, set())

    def test_sub_subset(self) -> None:
        nr1 = NumericRange(1, 15)
        nr2 = NumericRange(5, 10)
        nr = nr1 - nr2

        self.assertEqual(
            nr,
            {
                NumericRange(1, 4),
                NumericRange(11, 15),
            }
        )

    def test_xor_disjoint(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(15, 20)
        nr = nr1 ^ nr2

        self.assertEqual(
            nr,
            {
                NumericRange(1, 10),
                NumericRange(15, 20),
            }
        )

    def test_xor_adjacent(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(11, 20)
        nr = nr1 ^ nr2

        self.assertEqual(
            nr,
            {NumericRange(1, 20)}
        )

    def test_xor_equal(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(1, 10)
        nr = nr1 ^ nr2

        self.assertEqual(
            nr,
            set()
        )

    def test_xor(self) -> None:
        nr1 = NumericRange(1, 10)
        nr2 = NumericRange(5, 15)
        nr = nr1 ^ nr2

        self.assertEqual(
            nr,
            {
                NumericRange(1, 4),
                NumericRange(11, 15),
            }
        )

    def test_contains_outside(self) -> None:
        nr1 = NumericRange(1, 10)

        self.assertFalse(15 in nr1)

    def test_contains_outer_edge(self) -> None:
        nr1 = NumericRange(1, 10)

        self.assertFalse(0 in nr1)
        self.assertFalse(11 in nr1)

    def test_contains_inner_edge(self) -> None:
        nr1 = NumericRange(1, 10)

        self.assertTrue(1 in nr1)
        self.assertTrue(10 in nr1)

    def test_contains_inside(self) -> None:
        nr1 = NumericRange(1, 10)

        self.assertTrue(5 in nr1)

    def test_len(self) -> None:
        nr1 = NumericRange(1, 10)

        self.assertEqual(10, len(nr1))

    def test_len_minimum(self) -> None:
        nr1 = NumericRange(1, 1)

        self.assertEqual(1, len(nr1))


class TestNumericRangeSet(unittest.TestCase):

    def setUp(self) -> None:
        super().setUp()
        self.nrs = NumericRangeSet()

    def test_add_empty(self) -> None:
        with self.assertRaises(RuntimeError):
            self.nrs.add(10, 9)

    def test_add(self) -> None:
        self.nrs.add(1, 10)

        self.assertEqual(self.nrs.ranges, {NumericRange(1, 10)})

    def test_add_multiple_disjoint(self) -> None:
        self.nrs.add(1, 10)
        self.nrs.add(20, 30)

        self.assertEqual(
            self.nrs.ranges,
            {NumericRange(1, 10), NumericRange(20, 30)}
        )

    def test_add_multiple_adjacent(self) -> None:
        self.nrs.add(1, 10)
        self.nrs.add(11, 30)

        self.assertEqual(
            self.nrs.ranges,
            {NumericRange(1, 30)}
        )

    def test_add_multiple_overlap(self) -> None:
        self.nrs.add(1, 10)
        self.nrs.add(5, 30)

        self.assertEqual(
            self.nrs.ranges,
            {NumericRange(1, 30)}
        )

    def test_add_multiple_overlaps(self) -> None:
        self.nrs.add(5, 10)
        self.nrs.add(15, 20)
        self.nrs.add(25, 30)
        self.nrs.add(1, 30)

        self.assertEqual(
            self.nrs.ranges,
            {NumericRange(1, 30)}
        )

    def test_union_disjoint(self) -> None:
        self.nrs.add(15, 16)
        self.nrs.add(50, 60)

        nrs2 = NumericRangeSet()
        nrs2.add(1, 10)
        nrs2.add(20, 30)

        nrs_result = self.nrs | nrs2

        self.assertEqual(
            nrs_result.ranges,
            {
                NumericRange(1, 10),
                NumericRange(20, 30),
                NumericRange(15, 16),
                NumericRange(50, 60),
            }
        )

    def test_union_adjacent(self) -> None:
        self.nrs.add(15, 16)
        self.nrs.add(50, 60)

        nrs2 = NumericRangeSet()
        nrs2.add(1, 14)
        nrs2.add(20, 49)

        nrs_result = self.nrs | nrs2

        self.assertEqual(
            nrs_result.ranges,
            {
                NumericRange(1, 16),
                NumericRange(20, 60),
            }
        )

    def test_union_overlap(self) -> None:
        self.nrs.add(10, 16)
        self.nrs.add(50, 60)

        nrs2 = NumericRangeSet()
        nrs2.add(1, 14)
        nrs2.add(20, 60)

        nrs_result = self.nrs | nrs2

        self.assertEqual(
            nrs_result.ranges,
            {
                NumericRange(1, 16),
                NumericRange(20, 60),
            }
        )

    def test_intersection_empty_adjacent(self) -> None:
        self.nrs.add(15, 16)
        self.nrs.add(50, 60)

        nrs2 = NumericRangeSet()
        nrs2.add(1, 14)
        nrs2.add(20, 49)

        nrs_result = self.nrs & nrs2

        self.assertEqual(
            nrs_result.ranges,
            set()
        )

    def test_intersection_empty_disjoint(self) -> None:
        self.nrs.add(15, 16)
        self.nrs.add(50, 60)

        nrs2 = NumericRangeSet()
        nrs2.add(1, 10)
        nrs2.add(20, 40)

        nrs_result = self.nrs & nrs2

        self.assertEqual(
            nrs_result.ranges,
            set()
        )

    def test_intersection_overlap(self) -> None:
        self.nrs.add(1, 10)

        nrs2 = NumericRangeSet()
        nrs2.add(5, 15)

        nrs_result = self.nrs & nrs2

        self.assertEqual(
            nrs_result.ranges,
            {NumericRange(5, 10)}
        )

    def test_intersection_subset(self) -> None:
        self.nrs.add(1, 10)

        nrs2 = NumericRangeSet()
        nrs2.add(5, 9)

        nrs_result = self.nrs & nrs2

        self.assertEqual(
            nrs_result.ranges,
            {NumericRange(5, 9)}
        )

    def test_difference(self) -> None:
        self.nrs.add(1, 10)

        nrs2 = NumericRangeSet()
        nrs2.add(5, 15)

        nrs_result = self.nrs - nrs2

        self.assertEqual(
            nrs_result.ranges,
            {NumericRange(1, 4)}
        )

    def test_difference_disjoint(self) -> None:
        self.nrs.add(1, 10)

        nrs2 = NumericRangeSet()
        nrs2.add(13, 15)

        nrs_result = self.nrs - nrs2

        self.assertEqual(
            nrs_result.ranges,
            {NumericRange(1, 10)}
        )

    def test_difference_adjacent(self) -> None:
        self.nrs.add(1, 10)

        nrs2 = NumericRangeSet()
        nrs2.add(11, 15)

        nrs_result = self.nrs - nrs2

        self.assertEqual(
            nrs_result.ranges,
            {NumericRange(1, 10)}
        )

    def test_difference_subset_edge(self) -> None:
        self.nrs.add(1, 10)

        nrs2 = NumericRangeSet()
        nrs2.add(5, 10)

        nrs_result = self.nrs - nrs2

        self.assertEqual(
            nrs_result.ranges,
            {NumericRange(1, 4)}
        )

    def test_difference_subset_inner(self) -> None:
        self.nrs.add(1, 10)

        nrs2 = NumericRangeSet()
        nrs2.add(5, 7)

        nrs_result = self.nrs - nrs2

        self.assertEqual(
            nrs_result.ranges,
            {
                NumericRange(1, 4),
                NumericRange(8, 10),
            }
        )

    def test_symmetric_difference(self) -> None:
        self.nrs.add(1, 10)

        nrs2 = NumericRangeSet()
        nrs2.add(5, 15)

        nrs_result = self.nrs ^ nrs2

        self.assertEqual(
            nrs_result.ranges,
            {
                NumericRange(1, 4),
                NumericRange(11, 15),
            }
        )

    def test_symmetric_difference_disjoint(self) -> None:
        self.nrs.add(1, 10)

        nrs2 = NumericRangeSet()
        nrs2.add(12, 15)

        nrs_result = self.nrs ^ nrs2

        self.assertEqual(
            nrs_result.ranges,
            {
                NumericRange(1, 10),
                NumericRange(12, 15),
            }
        )

    def test_symmetric_difference_adjacent(self) -> None:
        self.nrs.add(1, 10)

        nrs2 = NumericRangeSet()
        nrs2.add(11, 15)

        nrs_result = self.nrs ^ nrs2

        self.assertEqual(
            nrs_result.ranges,
            {
                NumericRange(1, 15),
            }
        )

    def test_symmetric_difference_subset_edge(self) -> None:
        self.nrs.add(1, 10)

        nrs2 = NumericRangeSet()
        nrs2.add(7, 10)

        nrs_result = self.nrs ^ nrs2

        self.assertEqual(
            nrs_result.ranges,
            {
                NumericRange(1, 6),
            }
        )

    def test_symmetric_difference_subset_inner(self) -> None:
        self.nrs.add(1, 10)

        nrs2 = NumericRangeSet()
        nrs2.add(7, 8)

        nrs_result = self.nrs ^ nrs2

        self.assertEqual(
            nrs_result.ranges,
            {
                NumericRange(1, 6),
                NumericRange(9, 10),
            }
        )

    def test_contains_single_range(self) -> None:
        self.nrs.add(1, 10)

        self.assertTrue(1 in self.nrs)
        self.assertTrue(5 in self.nrs)
        self.assertTrue(10 in self.nrs)

    def test_not_contains_single_range(self) -> None:
        self.nrs.add(1, 10)

        self.assertFalse(0 in self.nrs)
        self.assertFalse(11 in self.nrs)
        self.assertFalse(15 in self.nrs)

    def test_contains_multiple_range(self) -> None:
        self.nrs.add(1, 10)
        self.nrs.add(20, 30)

        self.assertTrue(1 in self.nrs)
        self.assertTrue(5 in self.nrs)
        self.assertTrue(10 in self.nrs)
        self.assertTrue(20 in self.nrs)
        self.assertTrue(25 in self.nrs)
        self.assertTrue(30 in self.nrs)

    def test_not_contains_multiple_range(self) -> None:
        self.nrs.add(1, 10)
        self.nrs.add(20, 30)

        self.assertFalse(0 in self.nrs)
        self.assertFalse(11 in self.nrs)
        self.assertFalse(15 in self.nrs)
        self.assertFalse(19 in self.nrs)
        self.assertFalse(31 in self.nrs)
        self.assertFalse(35 in self.nrs)
