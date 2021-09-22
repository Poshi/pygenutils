'''pygenutils main module'''

import re
from dataclasses import dataclass
from functools import reduce, total_ordering
from typing import Dict, Iterator, Set

import pysam


class SequenceDict:

    def __init__(self, fasta: str):
        self.fasta = fasta
        self.sq: Dict[str, int] = dict()

        with pysam.FastaFile(fasta) as fin:  # pylint: disable=maybe-no-member
            self.sq = dict(zip(fin.references, fin.lengths))

    def __getitem__(self, chr) -> int:
        return self.sq[chr]

    def __contains__(self, chr) -> bool:
        return chr in self.sq

    def __iter__(self) -> Iterator:
        return iter(self.sq)

    def __len__(self) -> int:
        return len(self.sq)

    def __repr__(self) -> str:
        r = f"{self.__class__}({self.fasta})"

        return r

    def __str__(self) -> str:
        r = f"SequenceDict({self.fasta})"

        return r


@dataclass(frozen=True)
@total_ordering
class NumericRange:
    '''Helper class to represent a integer closed range of numbers.

    This class provides a basic imlementation of a single range of integers.
    It serves the purpose of doing the basic operations on ranges, considering
    them as sets. But it does not provide a uniform set of operations that work
    on this type of object and returns the same type of object.
    In particular, the only operation guaranteed to return a single range is the
    conjunction (__and__). The other operations can potentially return more than
    one single range, hence the result is a set of ranges instead of a single
    range.

    The class that provides full service for this kind of data is
    NumericRangeSet. This one works with sets of ranges and have the ability
    to work with any set of ranges and return any set of ranges.

    This class also provides some ordering between different ranges, so it
    contains most of the set operations but the __lt__ and others does not give
    information about subsets, but imposes an ordering between different
    instances of NumericRange.
    '''

    start: int
    end: int

    def __eq__(self, other) -> bool:
        if not isinstance(other, self.__class__):
            return NotImplemented

        return self.start == other.start and self.end == other.end

    def __lt__(self, other) -> bool:
        if not isinstance(other, self.__class__):
            return NotImplemented

        if self.start != other.start:
            return self.start < other.start

        return self.end < other.end

    def isdisjoint(self, other) -> bool:
        # If not the same class, cannot be compared
        if not isinstance(other, self.__class__):
            return NotImplemented

        # If it's the same class, we sort the tho objects according its start
        if self <= other:
            left, right = self, other
        else:
            left, right = other, self

        # The ranges are disjoint if left end comes before right start
        return left.end < right.start

    def isadjacent(self, other) -> bool:
        # If not the same class, cannot be compared
        if not isinstance(other, self.__class__):
            return NotImplemented

        # If it's the same class, we sort the tho objects according its start
        if self <= other:
            left, right = self, other
        else:
            left, right = other, self

        # The ranges are adjacent if left end comes right before right start
        return left.end + 1 == right.start

    def issubset(self, other) -> bool:
        return self.start >= other.start and self.end <= other.end

    def issuperset(self, other) -> bool:
        return self.start <= other.start and self.end >= other.end

    def __and__(self, other) -> 'NumericRange':
        if not isinstance(other, self.__class__):
            return NotImplemented

        if self.isdisjoint(other):
            # TODO: intersection of disjoint ranges should return an empty range, which cannot be represented here
            return NotImplemented

        return self.__class__(
            max(self.start, other.start),
            min(self.end, other.end),
        )

    def __or__(self, other) -> Set['NumericRange']:
        if not isinstance(other, self.__class__):
            return NotImplemented

        if self.isdisjoint(other) and not self.isadjacent(other):
            return {self, other}

        return {
            self.__class__(
                min(self.start, other.start),
                max(self.end, other.end),
            )
        }

    def __sub__(self, other) -> Set['NumericRange']:
        if not isinstance(other, self.__class__):
            return NotImplemented

        if self.isdisjoint(other):
            return {self}

        if other.issuperset(self):
            return set()

        if other.start <= self.start:
            return {NumericRange(other.end + 1, self.end)}
        elif other.end >= self.end:
            return {NumericRange(self.start, other.start - 1)}
        else:
            return {
                NumericRange(self.start, other.start - 1),
                NumericRange(other.end + 1, self.end),
            }

    def __xor__(self, other) -> Set['NumericRange']:
        if not isinstance(other, self.__class__):
            return NotImplemented

        if self.isadjacent(other):
            return self | other
        elif self.isdisjoint(other):
            return {self, other}
        elif self == other:
            return set()
        else:
            shared = self & other
            return {(self - shared).pop(), (other - shared).pop()}

    def __contains__(self, elem: int) -> bool:
        return self.start <= elem and elem <= self.end

    def __len__(self) -> int:
        return self.end - self.start + 1

    def __str__(self) -> str:
        r = f"[{self.start}, {self.end}]"

        return r

    def __repr__(self) -> str:
        r = f"{self.__class__}({self.start}, {self.end})"

        return r


class NumericRangeSet:

    def __init__(self, minimum: int, maximum: int):
        self.minimum: int = minimum
        self.maximum: int = maximum
        self.ranges: Set[NumericRange] = set()

    def add(self, start: int = None, end: int = None) -> 'NumericRangeSet':
        # Assign default values
        if start is None:
            start = self.minimum

        if end is None:
            end = self.maximum

        # Check values are inside the expected range
        if start < self.minimum or start >= self.maximum:
            raise RuntimeError(f'Start of range ({start}) not valid ([{self.minimum}, {self.maximum}])')

        if end < self.minimum or end >= self.maximum:
            raise RuntimeError(f'End of range ({end}) not valid ([{self.minimum}, {self.maximum}])')

        if start > end:
            raise RuntimeError(f'Start of range ({start}) greater than its end ({end})')

        # NumericRange is valid, so create it
        nr = NumericRange(start, end)

        # Extract the set of ranges that must be joined to the new NumericRange
        non_disjoint_ranges = {r for r in self.ranges if not r.isdisjoint(nr) or r.isadjacent(nr)}

        # Generate the new joined range and add it to the ranges
        def union_single_range(a: NumericRange, b: NumericRange) -> NumericRange:
            # We know for sure that the union will return a single item, so we
            # can just pop it from the set and return it for further processing.
            return (a | b).pop()
        self.ranges.add(reduce(union_single_range, non_disjoint_ranges, nr))

        # Remove the old ranges
        self.ranges -= non_disjoint_ranges

        return self

    def __and__(self, other) -> 'NumericRangeSet':
        # AKA conjunction
        if not isinstance(other, self.__class__):
            return NotImplemented

        if self.minimum != other.minimum or self.maximum != other.maximum:
            raise ValueError(
                f'Range boundaries does not match '
                f'({self.minimum}, {self.maximum})!=({other.minimum}, {other.maximum}). Cannot make union.'
            )

        result = NumericRangeSet(self.minimum, self.maximum)
        for r in self.ranges:
            intersect_ranges = [e for e in other.ranges if not r.isdisjoint(e)]
            result.ranges.update([e & r for e in intersect_ranges])

        return result

    def __or__(self, other) -> 'NumericRangeSet':
        # AKA union
        if not isinstance(other, self.__class__):
            return NotImplemented

        if self.minimum != other.minimum or self.maximum != other.maximum:
            raise ValueError(
                f'Range boundaries does not match '
                f'({self.minimum}, {self.maximum})!=({other.minimum}, {other.maximum}). Cannot make union.'
            )

        sorted_ranges = sorted(self.ranges | other.ranges, reverse=True)
        result_list = list()

        if sorted_ranges:
            result_list.append(sorted_ranges.pop())
        while sorted_ranges:
            last_range = result_list[-1]
            next_range = sorted_ranges.pop()
            if not next_range.isdisjoint(last_range) or next_range.isadjacent(last_range):
                result_list.pop()
                result_list.extend(next_range | last_range)
            else:
                result_list.append(next_range)

        result = NumericRangeSet(self.minimum, self.maximum)
        result.ranges = set(result_list)

        return result

    def __xor__(self, other) -> 'NumericRangeSet':
        # AKA symmetric difference
        return (self | other) - (self & other)

    def __sub__(self, other) -> 'NumericRangeSet':
        # AKA difference
        if not isinstance(other, self.__class__):
            return NotImplemented

        if self.minimum != other.minimum or self.maximum != other.maximum:
            raise ValueError(
                f'Range boundaries does not match '
                f'({self.minimum}, {self.maximum})!=({other.minimum}, {other.maximum}). Cannot make union.'
            )

        result = NumericRangeSet(self.minimum, self.maximum)
        for r in self.ranges:
            intersect_ranges = [e for e in other.ranges if not r.isdisjoint(e)]
            if not intersect_ranges:
                result.ranges.add(r)
            else:
                for e in intersect_ranges:
                    result.ranges.update(r - e)

        return result

    def __contains__(self, elem: int) -> bool:
        return any([elem in r for r in self.ranges])

    def __len__(self) -> int:
        return len(self.ranges)

    def __str__(self) -> str:
        r = f'NumericRangeSet({len(self.ranges)} ranges)'

        return r

    def __repr__(self) -> str:
        nrs = ", ".join([repr(nr) for nr in self.ranges])
        r = f'NumericRangeSet({self.minimum}, {self.maximum}, [{nrs}])'

        return r


class GenomicRangeSet:

    def __init__(self, sd: SequenceDict):
        self.sd = sd
        self.ranges: Dict[str, NumericRangeSet] = {k: NumericRangeSet(0, sd[k] - 1) for k in sd}
        self.string_re = re.compile(r'^(\w+)(:(\d+)?(-(\d+)?)?)?$')

    def add(self, chr: str, start: int = None, end: int = None) -> 'GenomicRangeSet':
        if chr not in self.ranges:
            raise ValueError(f'Chromosome {chr} not found in FastA file')

        self.ranges[chr].add(start, end)

        return self

    def add_from_string(self, s: str) -> 'GenomicRangeSet':
        '''Add interval from string.

        Parses the string and interprets it in the same way as samtools: one based index with both ends included.
        '''
        match = self.string_re.match(s)
        if match:
            chr, start, end = match.group(1, 3, 5)
            istart = int(start) - 1 if start is not None else None
            iend = int(end) - 1 if end is not None else None

            self.add(chr, istart, iend)
        else:
            raise ValueError(f'Bad formed region string: {s}')

        return self

    def __and__(self, other) -> 'GenomicRangeSet':
        if not isinstance(other, self.__class__):
            return NotImplemented

        if self.ranges.keys() != other.ranges.keys():
            raise ValueError("Keyset different in ranges, cannot make intersection")

        result = GenomicRangeSet(self.sd)
        for chr in result.ranges.keys():
            result.ranges[chr] = self.ranges[chr] & other.ranges[chr]

        return result

    def __or__(self, other) -> 'GenomicRangeSet':
        if not isinstance(other, self.__class__):
            return NotImplemented

        if self.ranges.keys() != other.ranges.keys():
            raise ValueError("Keyset different in ranges, cannot make union")

        result = GenomicRangeSet(self.sd)
        for chr in result.ranges.keys():
            result.ranges[chr] = self.ranges[chr] | other.ranges[chr]

        return result

    def __sub__(self, other) -> 'GenomicRangeSet':
        if not isinstance(other, self.__class__):
            return NotImplemented

        if self.ranges.keys() != other.ranges.keys():
            raise ValueError("Keyset different in ranges, cannot make union")

        result = GenomicRangeSet(self.sd)
        for chr in result.ranges.keys():
            result.ranges[chr] = self.ranges[chr] - other.ranges[chr]

        return result

    def __xor__(self, other) -> 'GenomicRangeSet':
        if not isinstance(other, self.__class__):
            return NotImplemented

        if self.ranges.keys() != other.ranges.keys():
            raise ValueError("Keyset different in ranges, cannot make union")

        result = GenomicRangeSet(self.sd)
        for chr in result.ranges.keys():
            result.ranges[chr] = self.ranges[chr] ^ other.ranges[chr]

        return result

    def __repr__(self) -> str:
        r = f'GenomicRangeSet({repr(self.ranges)})'

        return r

    def __str__(self) -> str:
        r = f'GenomicRangeSet({len(self.ranges)} chromosomes, {sum([len(v) for v in self.ranges.values()])} ranges)'

        return r


if __name__ == '__main__':
    print("Running in main")
    f = 'tests/data/test1.fasta'
    sd = SequenceDict(f)

    grs = GenomicRangeSet(sd)
