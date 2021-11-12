'''pygenutils main module'''

import locale
import re
from enum import Enum, auto, unique
from functools import reduce, total_ordering
from math import inf
from typing import Dict, Iterator, Literal, Optional, Set, Tuple, Union

from pysam import AlignedSegment  # pylint: disable=no-name-in-module
from pysam import AlignmentFile, FastaFile, VariantFile


class SequenceDict:
    '''Dictionary of sequences and their lengths.

    Represents a list of the names of the sequences and their associated
    lengths, so the ranges choosen can be validated.
    '''
    def __init__(self):
        self.seq_dict: Dict[str, int] = {}

    @classmethod
    def from_fasta(cls, fasta: str) -> 'SequenceDict':
        '''Create a SequenceDict from a FastA file'''
        result = cls()

        with FastaFile(fasta) as fin:
            result.seq_dict = dict(zip(fin.references, fin.lengths))

        return result

    @classmethod
    def from_bam(cls, bam: Union[str, AlignmentFile]) -> 'SequenceDict':
        '''Create a SequenceDict from a BAM file.

        The input can be either a string with the file name of a SAM or BAM file
        or an already opened AlignmentFile.
        '''
        result = cls()

        if isinstance(bam, AlignmentFile):
            sq_dict = bam.header['SQ']
        else:
            with AlignmentFile(bam) as fin:
                sq_dict = fin.header['SQ']

        result.seq_dict = {e['SN']: e['LN'] for e in sq_dict}

        return result

    @classmethod
    def from_bcf(cls, bcf: Union[str, VariantFile]) -> 'SequenceDict':
        '''Create a SequenceDict from a BCF file.

        The input can be either a string with the file name of a VCF or BCF file
        or an already opened VariantFile.
        '''
        result = cls()

        if isinstance(bcf, VariantFile):
            contigs = bcf.header.contigs
        else:
            with AlignmentFile(bcf) as fin:
                contigs = fin.header.contigs

        result.seq_dict = {e.name: e.length for e in contigs.values()}

        return result

    def __getitem__(self, chromosome) -> int:
        return self.seq_dict[chromosome]

    def __contains__(self, chromosome) -> bool:
        return chromosome in self.seq_dict

    def __iter__(self) -> Iterator[str]:
        return iter(self.seq_dict)

    def __len__(self) -> int:
        return len(self.seq_dict)

    def __repr__(self) -> str:
        return f"{self.__class__}({len(self.seq_dict)} sequences)"

    def __str__(self) -> str:
        return f"{self.__class__}({len(self.seq_dict)} sequences)"


NRElement = Union[int, float]


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

    def __init__(self, start: NRElement = -inf, end: NRElement = inf) -> None:
        if (
            (isinstance(start, float) and start != inf) or
            (isinstance(end, float) and end != inf)
        ):
            raise TypeError("No floats allowed in the range, only infinities")

        self.start = start
        self.end = end

    def as_tuple(self) -> Tuple[NRElement, NRElement]:
        '''Get the start and end elements in tuple format.'''
        return (self.start, self.end)

    def __hash__(self) -> int:
        return hash(self.as_tuple())

    def __eq__(self, other) -> bool:
        if not isinstance(other, self.__class__):
            return NotImplemented

        return self.start == other.start and self.end == other.end

    def __lt__(self, other) -> bool:
        if not isinstance(other, self.__class__):
            return NotImplemented

        if self.start == other.start:
            return self.end < other.end
        else:
            return self.start < other.start

    def isdisjoint(self, other) -> bool:
        '''Are the two intervals disjoint?

        Two intervals are disjoint when they share no positions. They can still
        be adjacent or not.
        '''
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
        '''Are the two intervals adjacent?

        Adjacent intervals are two intervals where the start position of one
        of them is innmediately after the end position of the other interval.
        '''
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
        '''Is one interval a subset of the other?'''
        # If not the same class, cannot be compared
        if not isinstance(other, self.__class__):
            return NotImplemented

        return self.start >= other.start and self.end <= other.end

    def issuperset(self, other) -> bool:
        '''Is one interval a superset of the other?'''
        # If not the same class, cannot be compared
        if not isinstance(other, self.__class__):
            return NotImplemented

        return self.start <= other.start and self.end >= other.end

    def __and__(self, other) -> Set['NumericRange']:
        if not isinstance(other, self.__class__):
            return NotImplemented

        if self.isdisjoint(other):
            return set()

        return {
            self.__class__(
                max(self.start, other.start),
                min(self.end, other.end),
            )
        }

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
            return {self.__class__(other.end + 1, self.end)}
        elif other.end >= self.end:
            return {self.__class__(self.start, other.start - 1)}
        else:
            return {
                self.__class__(self.start, other.start - 1),
                self.__class__(other.end + 1, self.end),
            }

    def __xor__(self, other) -> Set['NumericRange']:
        if not isinstance(other, self.__class__):
            return NotImplemented

        if self.isadjacent(other):
            return self | other

        if self.isdisjoint(other):
            return {self, other}

        if self == other:
            return set()

        # At this point we know that both intervals are overlapped, so the
        # shared region is a single one, and we can pop it.
        shared = (self & other).pop()

        return {(self - shared).pop(), (other - shared).pop()}

    def __contains__(self, elem: NRElement) -> bool:
        return self.start <= elem and elem <= self.end

    def __len__(self) -> int:
        # Float type is only used to store infinify. We check for the type so
        # the type checker knows that these variables are not integers
        if isinstance(self.start, float) or isinstance(self.end, float):
            raise ValueError('Infinite value length')

        return self.end - self.start + 1

    def __str__(self) -> str:
        return f"[{self.start}, {self.end}]"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.start}, {self.end})"


class NumericRangeSet:
    '''Represent a simple set of numeric intervals.

    This class allows us to perform set operations and existence testing on a
    set of numeric intervals, without any information about chromosomes.
    '''

    def __init__(self):
        self.ranges: Set[NumericRange] = set()

    def add(self, start: NRElement = -inf, end: NRElement = inf) -> 'NumericRangeSet':
        '''Add interval from its components'''
        if start is not None and end is not None and start > end:
            raise RuntimeError(f'Start of range ({start}) greater than its end ({end})')

        # Create the new object
        n_range = NumericRange(start, end)

        # Extract the set of ranges that must be joined to the new NumericRange
        non_disjoint_ranges = {
            r for r in self.ranges if not r.isdisjoint(n_range) or r.isadjacent(n_range)
        }

        # Generate the new joined range and add it to the ranges
        def union_single_range(nr1: NumericRange, nr2: NumericRange) -> NumericRange:
            # We know for sure that the union will return a single item, so we
            # can just pop it from the set and return it for further processing.
            return (nr1 | nr2).pop()
        self.ranges.add(reduce(union_single_range, non_disjoint_ranges, n_range))

        # Remove the old ranges
        self.ranges -= non_disjoint_ranges

        return self

    def __and__(self, other) -> 'NumericRangeSet':
        # AKA conjunction
        if not isinstance(other, self.__class__):
            return NotImplemented

        result = self.__class__()
        for rng in self.ranges:
            intersect_ranges = [e for e in other.ranges if not rng.isdisjoint(e)]
            result.ranges.update([(e & rng).pop() for e in intersect_ranges])

        return result

    def __or__(self, other) -> 'NumericRangeSet':
        # AKA union
        if not isinstance(other, self.__class__):
            return NotImplemented

        sorted_ranges = sorted(self.ranges | other.ranges, reverse=True)
        result_list = []

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

        result = self.__class__()
        result.ranges = set(result_list)

        return result

    def __xor__(self, other) -> 'NumericRangeSet':
        # AKA symmetric difference
        return (self | other) - (self & other)

    def __sub__(self, other) -> 'NumericRangeSet':
        # AKA difference
        if not isinstance(other, self.__class__):
            return NotImplemented

        result = self.__class__()
        for rng in self.ranges:
            intersect_ranges = [e for e in other.ranges if not rng.isdisjoint(e)]
            if not intersect_ranges:
                result.ranges.add(rng)
            else:
                for elem in intersect_ranges:
                    result.ranges.update(rng - elem)

        return result

    def __contains__(self, elem: NRElement) -> bool:
        return any(elem in r for r in self.ranges)

    def __iter__(self) -> Iterator[NumericRange]:
        return iter(self.ranges)

    def __len__(self) -> int:
        return len(self.ranges)

    def __str__(self) -> str:
        return f'{self.__class__.__name__}({len(self.ranges)} ranges)'

    def __repr__(self) -> str:
        nrs = ", ".join([repr(nr) for nr in self.ranges])
        return f'{self.__class__.__name__}([{nrs}])'


class GenomicRangeSet:
    '''Represents a fully operational set of intervals over a set of chromosomes.

    This class allows us to keep a list of chromosomes, each of the associated
    to a list of intervals. And be able to perform the whole set of set
    operations on them, together with existence testing.
    '''

    string_re = re.compile(r'^(\w+)(:(\d+)?(-(\d+)?)?)?$')

    def __init__(self, sequence_dict: SequenceDict = None):
        self.ranges: Dict[str, NumericRangeSet] = {}
        self.add_sequence_dict(sequence_dict)

    @staticmethod
    def _range_valid(
            sequence_dict: Optional[SequenceDict],
            chromosome: str,
            start: NRElement,
            end: NRElement
    ) -> None:
        '''Validate if a given range is inside the given SequenceDict.

        In case of an invalid range found, throw a ValueError exception,
        otherwise no operation is performed.
        If sequence_dict is None, all ranges are accepted
        '''

        if sequence_dict is not None:
            if chromosome not in sequence_dict:
                raise ValueError(f'Chromosome {chromosome} not found in the reference')
            if start != -inf and start < 0:
                raise ValueError(f'Start position below zero: {start}')
            if end != inf and end > sequence_dict[chromosome] - 1:
                raise ValueError(f'End position after end of reference: '
                                 f'{end + 1} > {sequence_dict[chromosome]}')

    def add_sequence_dict(self, sequence_dict: Optional[SequenceDict]) -> 'GenomicRangeSet':
        '''Add optional SequenceDict for validation of ranges.'''
        if sequence_dict is not None:
            # Validate current ranges, exception thrown if invalid
            for chromosome, nrs in self.ranges.items():
                for nrange in nrs:
                    self._range_valid(sequence_dict, chromosome, nrange.start, nrange.end)

        # Assign whatever we passed
        self.sequence_dict = sequence_dict

        return self

    def add(
        self,
        chromosome: str,
        start: NRElement = -inf,
        end: NRElement = inf
    ) -> 'GenomicRangeSet':
        '''Add interval from their components.'''
        # Validate data, exception thrown if invalid
        self._range_valid(self.sequence_dict, chromosome, start, end)

        # Data is OK, add the range
        if chromosome not in self.ranges:
            self.ranges[chromosome] = NumericRangeSet()

        self.ranges[chromosome].add(start, end)

        return self

    def add_from_string(self, s_range: str) -> 'GenomicRangeSet':
        '''Add interval from string.

        Parses the string and interprets it in the same way as samtools:
        one based index with both ends included.
        '''
        match = self.string_re.match(s_range)
        if match:
            chromosome, start, end = match.group(1, 3, 5)
            istart = int(start) - 1 if start is not None else -inf
            iend = int(end) - 1 if end is not None else inf

            self.add(chromosome, istart, iend)
        else:
            raise ValueError(f'Bad formed region string: {s_range}')

        return self

    @classmethod
    def from_bed_file(cls, bed_file_name: str, separator: str = None) -> 'GenomicRangeSet':
        '''Create a GenomicRangeSet that represents a BED file.'''
        result = cls()

        locale_encoding = locale.getpreferredencoding(False)
        with open(bed_file_name, 'r', encoding=locale_encoding) as bed_in:
            for line in bed_in:
                chromosome, start, end = line.rstrip().split(separator, maxsplit=3)[0:3]
                result.add(chromosome, int(start), int(end) - 1)

        return result

    def to_samtools_expression(self) -> str:
        '''Generate a string that contains a valid expression for filtering on a
        samtools view command.'''
        ranges_list = []
        for chromosome, nrs in self.ranges.items():
            for nrange in nrs:
                # BAM POS field is 1-based, so we have to add one to our numbers
                ranges_list.append(
                    f'(rname=~"{chromosome}" && pos>={nrange.start + 1} && pos<={nrange.end + 1})'
                )

        return ' || '.join(ranges_list)

    def to_bcftools_expression(self) -> str:
        '''Generate a string that contains a valid expression for filtering on a
        bcftools view command.'''
        ranges_list = []
        for chromosome, nrs in self.ranges.items():
            for nrange in nrs:
                # BCF POS field is 1-based, so we have to add one to our numbers
                ranges_list.append(
                    f'(CHROM=="{chromosome}" && POS>={nrange.start + 1} && POS<={nrange.end + 1})'
                )

        return ' || '.join(ranges_list)

    def to_regions(self) -> str:
        '''Generate a string that contains a valid range specification for the
        region parameter of samtools view and bcftools view.'''
        ranges_list = []
        for chromosome, nrs in self.ranges.items():
            for nrange in nrs:
                ranges_list.append(f'{chromosome}:{nrange.start + 1}-{nrange.end + 1}')

        return ' '.join(ranges_list)

    def __and__(self, other) -> 'GenomicRangeSet':
        if not isinstance(other, self.__class__):
            return NotImplemented

        result = self.__class__()
        for chromosome in self.ranges.keys() & other.ranges.keys():
            result.ranges[chromosome] = self.ranges[chromosome] & other.ranges[chromosome]

        return result

    def __or__(self, other) -> 'GenomicRangeSet':
        if not isinstance(other, self.__class__):
            return NotImplemented

        result = self.__class__()
        for chromosome in self.ranges.keys() & other.ranges.keys():
            result.ranges[chromosome] = self.ranges[chromosome] | other.ranges[chromosome]
        for chromosome in self.ranges.keys() - result.ranges.keys():
            result.ranges[chromosome] = self.ranges[chromosome]
        for chromosome in other.ranges.keys() - result.ranges.keys():
            result.ranges[chromosome] = other.ranges[chromosome]

        return result

    def __sub__(self, other) -> 'GenomicRangeSet':
        if not isinstance(other, self.__class__):
            return NotImplemented

        result = self.__class__()
        for chromosome in self.ranges.keys() & other.ranges.keys():
            result.ranges[chromosome] = self.ranges[chromosome] - other.ranges[chromosome]
        for chromosome in self.ranges.keys() - other.ranges.keys():
            result.ranges[chromosome] = self.ranges[chromosome]

        return result

    def __xor__(self, other) -> 'GenomicRangeSet':
        if not isinstance(other, self.__class__):
            return NotImplemented

        result = self.__class__()
        for chromosome in self.ranges.keys() & other.ranges.keys():
            result.ranges[chromosome] = self.ranges[chromosome] ^ other.ranges[chromosome]
        for chromosome in self.ranges.keys() - result.ranges.keys():
            result.ranges[chromosome] = self.ranges[chromosome]
        for chromosome in other.ranges.keys() - result.ranges.keys():
            result.ranges[chromosome] = other.ranges[chromosome]

        return result

    def __contains__(self, position: Tuple[str, NRElement]) -> bool:
        return position[0] in self.ranges and position[1] in self.ranges[position[0]]

    def __iter__(self) -> Iterator[str]:
        return iter(self.ranges)

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({repr(self.ranges)})'

    def __str__(self) -> str:
        return f'{self.__class__.__name__}({len(self.ranges)} chromosomes, ' \
               f'{sum([len(v) for v in self.ranges.values()])} ranges)'


@unique
class AlignmentFormat(Enum):
    '''Alignment file format identifiers.'''
    AUTO = auto()
    SAM = auto()
    BAM = auto()
    CRAM = auto()
    UNKNOWN = auto()

    @classmethod
    def guess_format(cls, filename, default=UNKNOWN) -> 'AlignmentFormat':
        '''Given a filename, guess in which format it is stored.'''
        if filename.endswith('.sam'):
            guessed = cls.SAM
        elif filename.endswith('.bam'):
            guessed = cls.BAM
        elif filename.endswith('.cram'):
            guessed = cls.CRAM
        else:
            guessed = default

        return guessed

    @classmethod
    def guess_mode_string(
            cls,
            filename: str,
            read_write: Literal['r', 'w'],
            default=UNKNOWN
    ) -> str:
        '''Utility method to compute the open file mode.'''
        file_format = cls.guess_format(filename, default)

        if file_format == cls.SAM:
            mode = ''
        elif file_format == cls.BAM:
            mode = 'b'
        elif file_format == cls.CRAM:
            mode = 'c'
        else:
            raise ValueError('Unknown alignment file format')

        return f'{read_write}{mode}'


@unique
class VariantFormat(Enum):
    '''Alignment file format identifiers.'''
    AUTO = auto()
    VCF = auto()
    BCF = auto()
    UNKNOWN = auto()

    @classmethod
    def guess_format(cls, filename, default=UNKNOWN) -> 'VariantFormat':
        '''Given a filename, guess in which format it is stored.'''
        if filename.endswith('.vcf'):
            guessed = cls.VCF
        elif filename.endswith('.vcf.gz'):
            guessed = cls.VCF
        elif filename.endswith('.bcf'):
            guessed = cls.BCF
        else:
            guessed = default

        return guessed

    @classmethod
    def guess_mode_string(
            cls,
            filename: str,
            read_write: Literal['r', 'w'],
            default=UNKNOWN
    ) -> str:
        '''Utility method to compute the open file mode.'''
        file_format = cls.guess_format(filename, default)

        if file_format == cls.VCF:
            mode = ''
        elif file_format == cls.BCF:
            mode = 'b'
        else:
            raise ValueError('Unknown variant file format')

        return f'{read_write}{mode}'


class BamFilter:
    '''Filter BAM records according their alignment positions.'''

    def __init__(
        self,
        grs: GenomicRangeSet,
        bam_file: AlignmentFile,
    ) -> None:
        self.grs = grs
        self.bam_file = bam_file

    def __iter__(self) -> AlignedSegment:
        for alignment in self.bam_file:
            # Reference start is 0-based leftmost coordinate
            if (alignment.reference_name, alignment.reference_start) in self.grs:
                yield alignment


class BcfFilter:
    '''Filter BCF records according their positions.'''

    def __init__(
        self,
        grs: GenomicRangeSet,
        bcf_file: VariantFile,
    ) -> None:
        self.grs = grs
        self.bcf_file = bcf_file

    def __iter__(self) -> AlignedSegment:
        for variant in self.bcf_file:
            # Position in variant files is 1-based and inclusive
            if (variant.chrom, variant.pos - 1) in self.grs:
                yield variant
