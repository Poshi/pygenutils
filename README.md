# pygenutils - Genomic Utilities for Python

This library pretends to be repository of utilities for working with genomics
data, covering for some holes found in the processing of data.
Many of these holes have been found when trying to process genomics data
in parallel by splitting the data in slices, with no alternatives found in
other commonly available tools.

To start with, it will include a couple clases to keep ranges of integers:

* NumericRangeSet
* GenomicRangeSet

NumericRangeSet intents to keep a set of ranges of integers, allowing set
operations between them.

GenomicRangeSet will encapsulate a set of NumericRangeSet together with
chromosome information, aiming to represent the intervals that could be
found in bed files.
