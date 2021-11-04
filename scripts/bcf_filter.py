#!python3

'''BCF filter.

Tool used to filter a BCF file and keep only all records whose starting position
is included in the specified range.
'''

import argparse
import sys
from contextlib import ExitStack

from pysam import VariantFile  # pylint: disable=no-name-in-module
from pygenutils import (VariantFormat, BcfFilter, GenomicRangeSet,
                        SequenceDict)


def process_clp() -> argparse.Namespace:
    '''Command line processing method'''
    parser = argparse.ArgumentParser(description='BAM filter')
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='BCF file to process',
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='resulting BCF file',
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

    return parser.parse_args()


if __name__ == '__main__':
    args = process_clp()

    # Check input parameters validity
    if not args.bed and not args.range:
        print('No valid ranges given', file=sys.stderr)
        sys.exit(1)

    # Open the contect manager to automatically close the opened files
    with ExitStack() as stack:
        # Open input file
        input_bcf = stack.enter_context(
            VariantFile(
                args.input,
                VariantFormat.guess_mode_string(args.input, 'r'),
            )
        )

        # Extract sequence dictionary for validation of the ranges
        sd = SequenceDict.from_bcf(input_bcf)

        # Load genomic ranges
        grs = GenomicRangeSet(sequence_dict=sd)
        for bed in args.bed:
            grs |= GenomicRangeSet.from_bed_file(bed)
        for rng in args.range:
            grs.add_from_string(rng)

        # Open output file
        output_bcf = stack.enter_context(
            VariantFile(
                args.output,
                VariantFormat.guess_mode_string(args.output, 'w'),
                header=input_bcf.header,
            )
        )

        # Finally, process the data
        bcf_filter = BcfFilter(grs, input_bcf)
        for variant in bcf_filter:
            output_bcf.write(variant)
