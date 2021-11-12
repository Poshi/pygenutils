#!python3

'''BAM filter.

Tool used to filter a BAM file and keep only all reads whose starting position
is included in the specified range.
'''

import argparse
import sys
from contextlib import ExitStack

from pysam import AlignmentFile  # pylint: disable=no-name-in-module

from pygenutils import (AlignmentFormat, BamFilter, GenomicRangeSet,
                        SequenceDict)


def process_clp() -> argparse.Namespace:
    '''Command line processing method'''
    parser = argparse.ArgumentParser(description='BAM filter')
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='BAM file to process',
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='resulting BAM file',
    )
    parser.add_argument(
        '-b', '--bed',
        required=False,
        action='append',
        default=[],
        help='BED file with ranges to include in output',
    )
    parser.add_argument(
        '-r', '--range',
        required=False,
        action='append',
        default=[],
        help='range to include in output',
    )
    parser.add_argument(
        '-R', '--reference',
        required=False,
        default=None,
        help='reference FastQ for CRAM operation',
    )

    return parser.parse_args()


def main():
    """Main bam_filter method"""
    args = process_clp()

    # Check input parameters validity
    if not args.bed and not args.range:
        print('No valid ranges given', file=sys.stderr)
        sys.exit(1)

    # Open the contect manager to automatically close the opened files
    with ExitStack() as stack:
        # Open input file
        input_bam = stack.enter_context(
            AlignmentFile(
                args.input,
                AlignmentFormat.guess_mode_string(args.input, 'r'),
                reference_filename=args.reference,
            )
        )

        # Extract sequence dictionary for validation of the ranges
        sequence_dict = SequenceDict.from_bam(input_bam)

        # Load genomic ranges
        grs = GenomicRangeSet(sequence_dict=sequence_dict)
        for bed in args.bed:
            grs |= GenomicRangeSet.from_bed_file(bed)
        for rng in args.range:
            grs.add_from_string(rng)

        # Open output file
        output_bam = stack.enter_context(
            AlignmentFile(
                args.output,
                AlignmentFormat.guess_mode_string(args.output, 'w'),
                template=input_bam,
                reference_filename=args.reference,
            )
        )

        # Finally, process the data
        bam_filter = BamFilter(grs, input_bam)
        for aln in bam_filter:
            output_bam.write(aln)


if __name__ == '__main__':
    sys.exit(main())
