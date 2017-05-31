

```bash
fasta_in="/input_data/mm9.fa"
gene_annot="/input_data/"
ann_dir="/input_data/"
ref_out="mm9"

```


```bash
#Bowtie
!bowtie
```

    bowtie2-build --threads 4 $fasta_in "$ann_dir""hg19"
    No output file specified!
    Bowtie 2 version 2.3.2 by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)
    Usage: bowtie2-build [options]* <reference_in> <bt2_index_base>
        reference_in            comma-separated list of files with ref sequences
        bt2_index_base          write bt2 data to files with this dir/basename
    *** Bowtie 2 indexes work only with v2 (not v1).  Likewise for v1 indexes. ***
    Options:
        -f                      reference files are Fasta (default)
        -c                      reference sequences given on cmd line (as
                                <reference_in>)
        --large-index           force generated index to be 'large', even if ref
                                has fewer than 4 billion nucleotides
        -a/--noauto             disable automatic -p/--bmax/--dcv memory-fitting
        -p/--packed             use packed strings internally; slower, less memory
        --bmax <int>            max bucket sz for blockwise suffix-array builder
        --bmaxdivn <int>        max bucket sz as divisor of ref len (default: 4)
        --dcv <int>             diff-cover period for blockwise (default: 1024)
        --nodc                  disable diff-cover (algorithm becomes quadratic)
        -r/--noref              don't build .3/.4 index files
        -3/--justref            just build .3/.4 index files
        -o/--offrate <int>      SA is sampled every 2^<int> BWT chars (default: 5)
        -t/--ftabchars <int>    # of chars consumed in initial lookup (default: 10)
        --threads <int>         # of threads
        --seed <int>            seed for random number generator
        -q/--quiet              verbose output (for debugging)
        -h/--help               print detailed description of tool and its options
        --usage                 print this usage message
        --version               print version information and quit





```bash
#Bowtie2 -- genome
bowtie2-build --threads 4 $fasta_in $ann_dir$ref_out
```

    Settings:
      Output files: "/input_data/mm9.*.bt2"
      Line rate: 6 (line is 64 bytes)
      Lines per side: 1 (side is 64 bytes)
      Offset rate: 4 (one in 16)
      FTable chars: 10
      Strings: unpacked
      Max bucket size: default
      Max bucket size, sqrt multiplier: default
      Max bucket size, len divisor: 16
      Difference-cover sample period: 1024
      Endianness: little
      Actual local endianness: little
      Sanity checking: disabled
      Assertions: disabled
      Random seed: 0
      Sizeofs: void*:8, int:4, long:8, size_t:8
    Input files DNA, FASTA:
      /input_data/mm9.fa
    Building a SMALL index
    Reading reference sizes
      Time reading reference sizes: 00:00:29
    Calculating joined length
    Writing header
    Reserving space for joined string
    Joining reference sequences
      Time to join reference sequences: 00:00:16
    bmax according to bmaxDivN setting: 163771623
    Using parameters --bmax 122828718 --dcv 1024
      Doing ahead-of-time memory usage test
      Passed!  Constructing with these parameters: --bmax 122828718 --dcv 1024
    Constructing suffix-array element generator
    Building DifferenceCoverSample
      Building sPrime
      Building sPrimeOrder
      V-Sorting samples
      V-Sorting samples time: 00:01:03
      Allocating rank array
      Ranking v-sort output
      Ranking v-sort output time: 00:00:37
      Invoking Larsson-Sadakane on ranks
      Invoking Larsson-Sadakane on ranks time: 00:00:59
      Sanity-checking and returning
    Building samples
    Reserving space for 44 sample suffixes
    Generating random suffixes
    QSorting 44 sample offsets, eliminating duplicates
    QSorting sample offsets, eliminating duplicates time: 00:00:00
    Multikey QSorting 44 samples
      (Using difference cover)
      Multikey QSorting samples time: 00:00:00
    Calculating bucket sizes
    Splitting and merging
      Splitting and merging time: 00:00:00
    Split 4, merged 21; iterating...
    Splitting and merging
      Splitting and merging time: 00:00:00
    Split 2, merged 2; iterating...
    Splitting and merging
      Splitting and merging time: 00:00:00
    Split 1, merged 1; iterating...
    Splitting and merging
      Splitting and merging time: 00:00:00
    Split 1, merged 0; iterating...
    Splitting and merging
      Splitting and merging time: 00:00:00
    Avg bucket size: 9.35838e+07 (target: 122828717)
    Converting suffix-array elements to index image
    Allocating ftab, absorbFtab
    Entering Ebwt loop
    Getting block 1 of 28
    Getting block 2 of 28
      Reserving size (122828718) for bucket 1
      Reserving size (122828718) for bucket 2
    Getting block 3 of 28
    Getting block 4 of 28
      Calculating Z arrays for bucket 1
      Calculating Z arrays for bucket 2
      Reserving size (122828718) for bucket 3
      Reserving size (122828718) for bucket 4
      Entering block accumulator loop for bucket 1:
      Entering block accumulator loop for bucket 2:
      Calculating Z arrays for bucket 3
      Calculating Z arrays for bucket 4
      Entering block accumulator loop for bucket 3:
      Entering block accumulator loop for bucket 4:
      bucket 1: 10%
      bucket 2: 10%
      bucket 3: 10%
      bucket 4: 10%
      bucket 1: 20%
      bucket 2: 20%
      bucket 3: 20%
      bucket 4: 20%
      bucket 1: 30%
      bucket 2: 30%
      bucket 3: 30%
      bucket 4: 30%
      bucket 1: 40%
      bucket 2: 40%
      bucket 1: 50%
      bucket 3: 40%
      bucket 4: 40%
      bucket 2: 50%
      bucket 1: 60%
      bucket 3: 50%
      bucket 4: 50%
      bucket 2: 60%
      bucket 1: 70%
      bucket 3: 60%
      bucket 4: 60%
      bucket 2: 70%
      bucket 1: 80%
      bucket 3: 70%
      bucket 2: 80%
      bucket 1: 90%
      bucket 4: 70%
      bucket 1: 100%
      Sorting block of length 69366514 for bucket 1
      (Using difference cover)
      bucket 2: 90%
      bucket 3: 80%
      bucket 4: 80%
      bucket 2: 100%
      Sorting block of length 107054429 for bucket 2
      (Using difference cover)
      bucket 3: 90%
      bucket 4: 90%
      bucket 3: 100%
      Sorting block of length 119483202 for bucket 3
      (Using difference cover)
      bucket 4: 100%
      Sorting block of length 76164057 for bucket 4
      (Using difference cover)
      Sorting block time: 00:01:21
    Returning block of 69366515 for bucket 1
    Getting block 5 of 28
      Reserving size (122828718) for bucket 5
      Calculating Z arrays for bucket 5
      Entering block accumulator loop for bucket 5:
      bucket 5: 10%
      Sorting block time: 00:01:34
    Returning block of 76164058 for bucket 4
      bucket 5: 20%
      bucket 5: 30%
      bucket 5: 40%
      bucket 5: 50%
      bucket 5: 60%
    Getting block 6 of 28
      Reserving size (122828718) for bucket 6
      Calculating Z arrays for bucket 6
      Entering block accumulator loop for bucket 6:
      bucket 5: 70%
      bucket 6: 10%
      bucket 5: 80%
      bucket 6: 20%
      bucket 5: 90%
      bucket 6: 30%
      bucket 5: 100%
      Sorting block of length 119927651 for bucket 5
      (Using difference cover)
      bucket 6: 40%
      bucket 6: 50%
      bucket 6: 60%
      Sorting block time: 00:02:14
    Returning block of 107054430 for bucket 2
      bucket 6: 70%
      bucket 6: 80%
      bucket 6: 90%
      bucket 6: 100%
      Sorting block of length 111329630 for bucket 6
      (Using difference cover)
      Sorting block time: 00:02:30
    Returning block of 119483203 for bucket 3
    Getting block 7 of 28
      Reserving size (122828718) for bucket 7
      Calculating Z arrays for bucket 7
      Entering block accumulator loop for bucket 7:
      bucket 7: 10%
      bucket 7: 20%
      bucket 7: 30%
      bucket 7: 40%
      bucket 7: 50%
      bucket 7: 60%
      bucket 7: 70%
    Getting block 8 of 28
      Reserving size (122828718) for bucket 8
      Calculating Z arrays for bucket 8
      Entering block accumulator loop for bucket 8:
      bucket 7: 80%
      bucket 8: 10%
      bucket 7: 90%
      bucket 8: 20%
      bucket 7: 100%
      Sorting block of length 56210234 for bucket 7
      (Using difference cover)
      bucket 8: 30%
      bucket 8: 40%
      bucket 8: 50%
      bucket 8: 60%
      bucket 8: 70%
      bucket 8: 80%
      bucket 8: 90%
      bucket 8: 100%
      Sorting block of length 71958774 for bucket 8
      (Using difference cover)
      Sorting block time: 00:01:24
    Returning block of 56210235 for bucket 7
    Getting block 9 of 28
      Reserving size (122828718) for bucket 9
      Calculating Z arrays for bucket 9
      Entering block accumulator loop for bucket 9:
      bucket 9: 10%
      bucket 9: 20%
      bucket 9: 30%
      Sorting block time: 00:02:56
    Returning block of 119927652 for bucket 5
      bucket 9: 40%
      bucket 9: 50%
      Sorting block time: 00:02:42
    Returning block of 111329631 for bucket 6
      bucket 9: 60%
      bucket 9: 70%
      bucket 9: 80%
      Sorting block time: 00:01:48
    Returning block of 71958775 for bucket 8
      bucket 9: 90%
      bucket 9: 100%
      Sorting block of length 116909663 for bucket 9
      (Using difference cover)
    Getting block 10 of 28
      Reserving size (122828718) for bucket 10
      Calculating Z arrays for bucket 10
      Entering block accumulator loop for bucket 10:
    Getting block 11 of 28
      Reserving size (122828718) for bucket 11
      Calculating Z arrays for bucket 11
      Entering block accumulator loop for bucket 11:
      bucket 10: 10%
      bucket 11: 10%
    Getting block 12 of 28
      Reserving size (122828718) for bucket 12
      Calculating Z arrays for bucket 12
      Entering block accumulator loop for bucket 12:
      bucket 10: 20%
      bucket 11: 20%
      bucket 12: 10%
      bucket 10: 30%
      bucket 11: 30%
      bucket 12: 20%
      bucket 10: 40%
      bucket 11: 40%
      bucket 12: 30%
      bucket 10: 50%
      bucket 11: 50%
      bucket 12: 40%
      bucket 10: 60%
      bucket 11: 60%
      bucket 10: 70%
      bucket 12: 50%
      bucket 11: 70%
      bucket 10: 80%
      bucket 12: 60%
      bucket 11: 80%
      bucket 10: 90%
      bucket 12: 70%
      bucket 11: 90%
      bucket 10: 100%
      Sorting block of length 71463388 for bucket 10
      (Using difference cover)
      bucket 12: 80%
      bucket 11: 100%
      Sorting block of length 122449757 for bucket 11
      (Using difference cover)
      bucket 12: 90%
      bucket 12: 100%
      Sorting block of length 106198633 for bucket 12
      (Using difference cover)
      Sorting block time: 00:01:34
    Returning block of 71463389 for bucket 10
    Getting block 13 of 28
      Reserving size (122828718) for bucket 13
      Calculating Z arrays for bucket 13
      Entering block accumulator loop for bucket 13:
      bucket 13: 10%
      bucket 13: 20%
      Sorting block time: 00:02:40
    Returning block of 116909664 for bucket 9
      bucket 13: 30%
      bucket 13: 40%
      bucket 13: 50%
      bucket 13: 60%
      bucket 13: 70%
      bucket 13: 80%
      bucket 13: 90%
    Getting block 14 of 28
      Reserving size (122828718) for bucket 14
      Calculating Z arrays for bucket 14
      Entering block accumulator loop for bucket 14:
      bucket 13: 100%
      Sorting block of length 114547269 for bucket 13
      (Using difference cover)
      bucket 14: 10%
      bucket 14: 20%
      Sorting block time: 00:02:24
    Returning block of 106198634 for bucket 12
      bucket 14: 30%
      bucket 14: 40%
      bucket 14: 50%
      bucket 14: 60%
      Sorting block time: 00:02:48
    Returning block of 122449758 for bucket 11
      bucket 14: 70%
      bucket 14: 80%
    Getting block 15 of 28
      Reserving size (122828718) for bucket 15
      Calculating Z arrays for bucket 15
      Entering block accumulator loop for bucket 15:
      bucket 14: 90%
      bucket 15: 10%
      bucket 14: 100%
      Sorting block of length 101495912 for bucket 14
      (Using difference cover)
      bucket 15: 20%
      bucket 15: 30%
      bucket 15: 40%
      bucket 15: 50%
    Getting block 16 of 28
      Reserving size (122828718) for bucket 16
      Calculating Z arrays for bucket 16
      Entering block accumulator loop for bucket 16:
      bucket 15: 60%
      bucket 16: 10%
      bucket 15: 70%
      bucket 16: 20%
      bucket 15: 80%
      bucket 16: 30%
      bucket 15: 90%
      bucket 16: 40%
      bucket 15: 100%
      Sorting block of length 72277302 for bucket 15
      (Using difference cover)
      bucket 16: 50%
      bucket 16: 60%
      bucket 16: 70%
      bucket 16: 80%
      bucket 16: 90%
      bucket 16: 100%
      Sorting block of length 112423886 for bucket 16
      (Using difference cover)
      Sorting block time: 00:02:36
    Returning block of 114547270 for bucket 13
      Sorting block time: 00:01:32
    Returning block of 72277303 for bucket 15
      Sorting block time: 00:02:19
    Returning block of 101495913 for bucket 14
    Getting block 17 of 28
      Reserving size (122828718) for bucket 17
      Calculating Z arrays for bucket 17
      Entering block accumulator loop for bucket 17:
    Getting block 18 of 28
      Reserving size (122828718) for bucket 18
      Calculating Z arrays for bucket 18
      Entering block accumulator loop for bucket 18:
      bucket 17: 10%
      bucket 18: 10%
      bucket 18: 20%
      bucket 17: 20%
      bucket 18: 30%
      bucket 17: 30%
      bucket 18: 40%
      bucket 17: 40%
      bucket 18: 50%
      bucket 17: 50%
    Getting block 19 of 28
      Reserving size (122828718) for bucket 19
      Calculating Z arrays for bucket 19
      Entering block accumulator loop for bucket 19:
      bucket 18: 60%
      bucket 17: 60%
      bucket 19: 10%
      bucket 18: 70%
      bucket 17: 70%
      bucket 19: 20%
      bucket 18: 80%
      bucket 17: 80%
      bucket 19: 30%
      bucket 18: 90%
      bucket 17: 90%
      bucket 19: 40%
      bucket 18: 100%
      Sorting block of length 100679308 for bucket 18
      (Using difference cover)
      bucket 17: 100%
      Sorting block of length 113135159 for bucket 17
      (Using difference cover)
      bucket 19: 50%
      bucket 19: 60%
      bucket 19: 70%
      bucket 19: 80%
      bucket 19: 90%
      bucket 19: 100%
      Sorting block of length 57786489 for bucket 19
      (Using difference cover)
      Sorting block time: 00:02:36
    Returning block of 112423887 for bucket 16
    Getting block 20 of 28
      Reserving size (122828718) for bucket 20
      Calculating Z arrays for bucket 20
      Entering block accumulator loop for bucket 20:
      bucket 20: 10%
      bucket 20: 20%
      bucket 20: 30%
      bucket 20: 40%
      bucket 20: 50%
      bucket 20: 60%
      bucket 20: 70%
      bucket 20: 80%
      bucket 20: 90%
      bucket 20: 100%
      Sorting block of length 96659056 for bucket 20
      (Using difference cover)
      Sorting block time: 00:01:21
    Returning block of 57786490 for bucket 19
    Getting block 21 of 28
      Reserving size (122828718) for bucket 21
      Calculating Z arrays for bucket 21
      Entering block accumulator loop for bucket 21:
      bucket 21: 10%
      bucket 21: 20%
      bucket 21: 30%
      bucket 21: 40%
      bucket 21: 50%
      bucket 21: 60%
      bucket 21: 70%
      bucket 21: 80%
      bucket 21: 90%
      Sorting block time: 00:02:30
    Returning block of 100679309 for bucket 18
      bucket 21: 100%
      Sorting block of length 29805507 for bucket 21
      (Using difference cover)
      Sorting block time: 00:02:55
    Returning block of 113135160 for bucket 17
    Getting block 22 of 28
      Reserving size (122828718) for bucket 22
      Calculating Z arrays for bucket 22
      Entering block accumulator loop for bucket 22:
      bucket 22: 10%
      bucket 22: 20%
      bucket 22: 30%
      bucket 22: 40%
      bucket 22: 50%
      bucket 22: 60%
      Sorting block time: 00:00:49
    Returning block of 29805508 for bucket 21
    Getting block 23 of 28
      Reserving size (122828718) for bucket 23
      Calculating Z arrays for bucket 23
      Entering block accumulator loop for bucket 23:
      bucket 22: 70%
      bucket 23: 10%
      bucket 22: 80%
    Getting block 24 of 28
      Reserving size (122828718) for bucket 24
      Calculating Z arrays for bucket 24
      Entering block accumulator loop for bucket 24:
      bucket 23: 20%
      bucket 22: 90%
      bucket 24: 10%
      bucket 23: 30%
      bucket 22: 100%
      Sorting block of length 111607091 for bucket 22
      (Using difference cover)
      bucket 24: 20%
      bucket 23: 40%
      bucket 24: 30%
      bucket 23: 50%
      bucket 24: 40%
      bucket 23: 60%
      bucket 24: 50%
      bucket 23: 70%
      bucket 24: 60%
      bucket 23: 80%
      bucket 24: 70%
      bucket 23: 90%
      bucket 23: 100%
      Sorting block of length 44391035 for bucket 23
      (Using difference cover)
      bucket 24: 80%
      bucket 24: 90%
      Sorting block time: 00:02:29
    Returning block of 96659057 for bucket 20
      bucket 24: 100%
      Sorting block of length 106990899 for bucket 24
      (Using difference cover)
    Getting block 25 of 28
      Reserving size (122828718) for bucket 25
      Calculating Z arrays for bucket 25
      Entering block accumulator loop for bucket 25:
      bucket 25: 10%
      bucket 25: 20%
      bucket 25: 30%
      bucket 25: 40%
      bucket 25: 50%
      bucket 25: 60%
      bucket 25: 70%
      bucket 25: 80%
      Sorting block time: 00:01:13
    Returning block of 44391036 for bucket 23
      bucket 25: 90%
      bucket 25: 100%
      Sorting block of length 95171746 for bucket 25
      (Using difference cover)
    Getting block 26 of 28
      Reserving size (122828718) for bucket 26
      Calculating Z arrays for bucket 26
      Entering block accumulator loop for bucket 26:
      bucket 26: 10%
      bucket 26: 20%
      bucket 26: 30%
      bucket 26: 40%
      bucket 26: 50%
      bucket 26: 60%
      bucket 26: 70%
      bucket 26: 80%
      bucket 26: 90%
      bucket 26: 100%
      Sorting block of length 95802932 for bucket 26
      (Using difference cover)
      Sorting block time: 00:02:50
    Returning block of 111607092 for bucket 22
    Getting block 27 of 28
      Reserving size (122828718) for bucket 27
      Calculating Z arrays for bucket 27
      Entering block accumulator loop for bucket 27:
      bucket 27: 10%
      bucket 27: 20%
      Sorting block time: 00:03:01
    Returning block of 106990900 for bucket 24
      bucket 27: 30%
      bucket 27: 40%
      bucket 27: 50%
      bucket 27: 60%
      bucket 27: 70%
      bucket 27: 80%
    Getting block 28 of 28
      Reserving size (122828718) for bucket 28
      Calculating Z arrays for bucket 28
      Entering block accumulator loop for bucket 28:
      bucket 27: 90%
      bucket 28: 10%
      bucket 27: 100%
      Sorting block of length 120049009 for bucket 27
      (Using difference cover)
      bucket 28: 20%
      bucket 28: 30%
      bucket 28: 40%
      bucket 28: 50%
      bucket 28: 60%
      bucket 28: 70%
      bucket 28: 80%
      bucket 28: 90%
      Sorting block time: 00:02:53
    Returning block of 95171747 for bucket 25
      bucket 28: 100%
      Sorting block of length 99007413 for bucket 28
      (Using difference cover)
      Sorting block time: 00:03:07
    Returning block of 95802933 for bucket 26
      Sorting block time: 00:02:46
    Returning block of 99007414 for bucket 28
      Sorting block time: 00:03:26
    Returning block of 120049010 for bucket 27
    Exited Ebwt loop
    fchr[A]: 0
    fchr[C]: 763450008
    fchr[G]: 1309739408
    fchr[T]: 1856187599
    fchr[$]: 2620345972
    Exiting Ebwt::buildToDisk()
    Returning from initFromVector
    Wrote 877656917 bytes to primary EBWT file: /input_data/mm9.1.bt2
    Wrote 655086500 bytes to secondary EBWT file: /input_data/mm9.2.bt2
    Re-opening _in1 and _in2 as input streams
    Returning from Ebwt constructor
    Headers:
        len: 2620345972
        bwtLen: 2620345973
        sz: 655086493
        bwtSz: 655086494
        lineRate: 6
        offRate: 4
        offMask: 0xfffffff0
        ftabChars: 10
        eftabLen: 20
        eftabSz: 80
        ftabLen: 1048577
        ftabSz: 4194308
        offsLen: 163771624
        offsSz: 655086496
        lineSz: 64
        sideSz: 64
        sideBwtSz: 48
        sideBwtLen: 192
        numSides: 13647636
        numLines: 13647636
        ebwtTotLen: 873448704
        ebwtTotSz: 873448704
        color: 0
        reverse: 0
    Total time for call to driver() for forward index: 00:31:50
    Reading reference sizes
      Time reading reference sizes: 00:00:32
    Calculating joined length
    Writing header
    Reserving space for joined string
    Joining reference sequences
      Time to join reference sequences: 00:00:15
      Time to reverse reference sequence: 00:00:02
    bmax according to bmaxDivN setting: 163771623
    Using parameters --bmax 122828718 --dcv 1024
      Doing ahead-of-time memory usage test
      Passed!  Constructing with these parameters: --bmax 122828718 --dcv 1024
    Constructing suffix-array element generator
    Building DifferenceCoverSample
      Building sPrime
      Building sPrimeOrder
      V-Sorting samples
      V-Sorting samples time: 00:01:13
      Allocating rank array
      Ranking v-sort output
      Ranking v-sort output time: 00:00:36
      Invoking Larsson-Sadakane on ranks
      Invoking Larsson-Sadakane on ranks time: 00:01:03
      Sanity-checking and returning
    Building samples
    Reserving space for 44 sample suffixes
    Generating random suffixes
    QSorting 44 sample offsets, eliminating duplicates
    QSorting sample offsets, eliminating duplicates time: 00:00:00
    Multikey QSorting 44 samples
      (Using difference cover)
      Multikey QSorting samples time: 00:00:00
    Calculating bucket sizes
    Splitting and merging
      Splitting and merging time: 00:00:00
    Split 4, merged 20; iterating...
    Splitting and merging
      Splitting and merging time: 00:00:00
    Split 3, merged 2; iterating...
    Splitting and merging
      Splitting and merging time: 00:00:00
    Split 2, merged 2; iterating...
    Splitting and merging
      Splitting and merging time: 00:00:00
    Avg bucket size: 9.35838e+07 (target: 122828717)
    Converting suffix-array elements to index image
    Allocating ftab, absorbFtab
    Entering Ebwt loop
    Getting block 1 of 28
      Reserving size (122828718) for bucket 1
    Getting block 3 of 28
    Getting block 2 of 28
      Calculating Z arrays for bucket 1
      Reserving size (122828718) for bucket 3
      Reserving size (122828718) for bucket 2
      Calculating Z arrays for bucket 3
      Calculating Z arrays for bucket 2
      Entering block accumulator loop for bucket 3:
      Entering block accumulator loop for bucket 1:
      Entering block accumulator loop for bucket 2:
    Getting block 4 of 28
      Reserving size (122828718) for bucket 4
      Calculating Z arrays for bucket 4
      Entering block accumulator loop for bucket 4:
      bucket 1: 10%
      bucket 2: 10%
      bucket 3: 10%
      bucket 4: 10%
      bucket 1: 20%
      bucket 2: 20%
      bucket 3: 20%
      bucket 4: 20%
      bucket 1: 30%
      bucket 2: 30%
      bucket 3: 30%
      bucket 4: 30%
      bucket 1: 40%
      bucket 2: 40%
      bucket 3: 40%
      bucket 4: 40%
      bucket 1: 50%
      bucket 2: 50%
      bucket 3: 50%
      bucket 1: 60%
      bucket 4: 50%
      bucket 2: 60%
      bucket 1: 70%
      bucket 3: 60%
      bucket 4: 60%
      bucket 2: 70%
      bucket 1: 80%
      bucket 3: 70%
      bucket 2: 80%
      bucket 4: 70%
      bucket 1: 90%
      bucket 2: 90%
      bucket 3: 80%
      bucket 4: 80%
      bucket 1: 100%
      Sorting block of length 117942128 for bucket 1
      (Using difference cover)
      bucket 2: 100%
      Sorting block of length 103666198 for bucket 2
      (Using difference cover)
      bucket 3: 90%
      bucket 4: 90%
      bucket 3: 100%
      Sorting block of length 72732200 for bucket 3
      (Using difference cover)
      bucket 4: 100%
      Sorting block of length 66774133 for bucket 4
      (Using difference cover)
      Sorting block time: 00:01:25
    Returning block of 66774134 for bucket 4
      Sorting block time: 00:01:29
    Returning block of 72732201 for bucket 3
    Getting block 5 of 28
      Reserving size (122828718) for bucket 5
      Calculating Z arrays for bucket 5
      Entering block accumulator loop for bucket 5:
      bucket 5: 10%
    Getting block 6 of 28
      Reserving size (122828718) for bucket 6
      Calculating Z arrays for bucket 6
      Entering block accumulator loop for bucket 6:
      bucket 5: 20%
      bucket 6: 10%
      bucket 5: 30%
      bucket 6: 20%
      bucket 5: 40%
      bucket 6: 30%
      bucket 5: 50%
      bucket 6: 40%
      bucket 5: 60%
      bucket 6: 50%
      bucket 5: 70%
      bucket 6: 60%
      bucket 5: 80%
      bucket 6: 70%
      bucket 5: 90%
      bucket 6: 80%
      bucket 5: 100%
      Sorting block of length 114557282 for bucket 5
      (Using difference cover)
      bucket 6: 90%
      bucket 6: 100%
      Sorting block of length 115432025 for bucket 6
      (Using difference cover)
      Sorting block time: 00:02:21
    Returning block of 103666199 for bucket 2
      Sorting block time: 00:02:40
    Returning block of 117942129 for bucket 1
    Getting block 7 of 28
      Reserving size (122828718) for bucket 7
      Calculating Z arrays for bucket 7
      Entering block accumulator loop for bucket 7:
      bucket 7: 10%
      bucket 7: 20%
      bucket 7: 30%
      bucket 7: 40%
      bucket 7: 50%
      bucket 7: 60%
    Getting block 8 of 28
      Reserving size (122828718) for bucket 8
      Calculating Z arrays for bucket 8
      Entering block accumulator loop for bucket 8:
      bucket 7: 70%
      bucket 8: 10%
      bucket 7: 80%
      bucket 8: 20%
      bucket 7: 90%
      bucket 8: 30%
      bucket 7: 100%
      Sorting block of length 98740313 for bucket 7
      (Using difference cover)
      bucket 8: 40%
      bucket 8: 50%
      bucket 8: 60%
      bucket 8: 70%
      bucket 8: 80%
      bucket 8: 90%
      bucket 8: 100%
      Sorting block of length 120687041 for bucket 8
      (Using difference cover)
      Sorting block time: 00:02:41
    Returning block of 114557283 for bucket 5
      Sorting block time: 00:02:44
    Returning block of 115432026 for bucket 6
    Getting block 9 of 28
      Reserving size (122828718) for bucket 9
      Calculating Z arrays for bucket 9
      Entering block accumulator loop for bucket 9:
      bucket 9: 10%
      bucket 9: 20%
    Getting block 10 of 28
      Reserving size (122828718) for bucket 10
      Calculating Z arrays for bucket 10
      Entering block accumulator loop for bucket 10:
      bucket 9: 30%
      bucket 10: 10%
      Sorting block time: 00:02:18
    Returning block of 98740314 for bucket 7
      bucket 9: 40%
      bucket 10: 20%
      bucket 9: 50%
      bucket 10: 30%
      bucket 9: 60%
      bucket 10: 40%
      bucket 9: 70%
      bucket 10: 50%
      bucket 9: 80%
      bucket 10: 60%
      bucket 9: 90%
      bucket 10: 70%
      bucket 9: 100%
      Sorting block of length 98799048 for bucket 9
      (Using difference cover)
    Getting block 11 of 28
      Reserving size (122828718) for bucket 11
      Calculating Z arrays for bucket 11
      Entering block accumulator loop for bucket 11:
      bucket 10: 80%
      bucket 10: 90%
      bucket 11: 10%
      bucket 10: 100%
      Sorting block of length 106167244 for bucket 10
      (Using difference cover)
      bucket 11: 20%
      bucket 11: 30%
      bucket 11: 40%
      bucket 11: 50%
      bucket 11: 60%
      bucket 11: 70%
      bucket 11: 80%
      bucket 11: 90%
      Sorting block time: 00:02:53
    Returning block of 120687042 for bucket 8
      bucket 11: 100%
      Sorting block of length 60764740 for bucket 11
      (Using difference cover)
    Getting block 12 of 28
      Reserving size (122828718) for bucket 12
      Calculating Z arrays for bucket 12
      Entering block accumulator loop for bucket 12:
      bucket 12: 10%
      bucket 12: 20%
      bucket 12: 30%
      bucket 12: 40%
      bucket 12: 50%
      bucket 12: 60%
      bucket 12: 70%
      bucket 12: 80%
      bucket 12: 90%
      bucket 12: 100%
      Sorting block of length 104307603 for bucket 12
      (Using difference cover)
      Sorting block time: 00:01:23
    Returning block of 60764741 for bucket 11
    Getting block 13 of 28
      Reserving size (122828718) for bucket 13
      Calculating Z arrays for bucket 13
      Entering block accumulator loop for bucket 13:
      Sorting block time: 00:02:17
    Returning block of 98799049 for bucket 9
      bucket 13: 10%
      bucket 13: 20%
      bucket 13: 30%
      bucket 13: 40%
      bucket 13: 50%
    Getting block 14 of 28
      Reserving size (122828718) for bucket 14
      Calculating Z arrays for bucket 14
      Entering block accumulator loop for bucket 14:
      bucket 13: 60%
      Sorting block time: 00:02:30
    Returning block of 106167245 for bucket 10
      bucket 14: 10%
      bucket 13: 70%
      bucket 14: 20%
      bucket 13: 80%
      bucket 14: 30%
      bucket 13: 90%
      bucket 14: 40%
      bucket 13: 100%
      Sorting block of length 118148363 for bucket 13
      (Using difference cover)
      bucket 14: 50%
      bucket 14: 60%
    Getting block 15 of 28
      Reserving size (122828718) for bucket 15
      Calculating Z arrays for bucket 15
      Entering block accumulator loop for bucket 15:
      bucket 14: 70%
      bucket 15: 10%
      bucket 14: 80%
      bucket 15: 20%
      bucket 14: 90%
      bucket 15: 30%
      bucket 14: 100%
      Sorting block of length 72526733 for bucket 14
      (Using difference cover)
      bucket 15: 40%
      bucket 15: 50%
      bucket 15: 60%
      bucket 15: 70%
      bucket 15: 80%
      bucket 15: 90%
      bucket 15: 100%
      Sorting block of length 51302171 for bucket 15
      (Using difference cover)
      Sorting block time: 00:02:36
    Returning block of 104307604 for bucket 12
    Getting block 16 of 28
      Reserving size (122828718) for bucket 16
      Calculating Z arrays for bucket 16
      Entering block accumulator loop for bucket 16:
      bucket 16: 10%
      bucket 16: 20%
      bucket 16: 30%
      bucket 16: 40%
      Sorting block time: 00:01:21
    Returning block of 51302172 for bucket 15
      bucket 16: 50%
      bucket 16: 60%
      bucket 16: 70%
      Sorting block time: 00:01:56
    Returning block of 72526734 for bucket 14
      bucket 16: 80%
    Getting block 17 of 28
      Reserving size (122828718) for bucket 17
      Calculating Z arrays for bucket 17
      Entering block accumulator loop for bucket 17:
      bucket 16: 90%
      bucket 17: 10%
      bucket 16: 100%
      Sorting block of length 111550527 for bucket 16
      (Using difference cover)
      bucket 17: 20%
      bucket 17: 30%
      bucket 17: 40%
    Getting block 18 of 28
      Reserving size (122828718) for bucket 18
      Calculating Z arrays for bucket 18
      Entering block accumulator loop for bucket 18:
      bucket 17: 50%
      bucket 18: 10%
      bucket 17: 60%
      bucket 18: 20%
      bucket 17: 70%
      bucket 18: 30%
      bucket 17: 80%
      bucket 18: 40%
      bucket 17: 90%
      bucket 18: 50%
      bucket 17: 100%
      Sorting block of length 100808090 for bucket 17
      (Using difference cover)
      bucket 18: 60%
      bucket 18: 70%
      Sorting block time: 00:03:06
    Returning block of 118148364 for bucket 13
      bucket 18: 80%
      bucket 18: 90%
      bucket 18: 100%
      Sorting block of length 80329239 for bucket 18
      (Using difference cover)
    Getting block 19 of 28
      Reserving size (122828718) for bucket 19
      Calculating Z arrays for bucket 19
      Entering block accumulator loop for bucket 19:
      bucket 19: 10%
      bucket 19: 20%
      bucket 19: 30%
      bucket 19: 40%
      bucket 19: 50%
      bucket 19: 60%
      bucket 19: 70%
      bucket 19: 80%
      bucket 19: 90%
      bucket 19: 100%
      Sorting block of length 47393288 for bucket 19
      (Using difference cover)
      Sorting block time: 00:02:09
    Returning block of 80329240 for bucket 18
      Sorting block time: 00:02:59
    Returning block of 111550528 for bucket 16
      Sorting block time: 00:01:17
    Returning block of 47393289 for bucket 19
      Sorting block time: 00:02:49
    Returning block of 100808091 for bucket 17
    Getting block 20 of 28
      Reserving size (122828718) for bucket 20
      Calculating Z arrays for bucket 20
      Entering block accumulator loop for bucket 20:
      bucket 20: 10%
    Getting block 21 of 28
      Reserving size (122828718) for bucket 21
      Calculating Z arrays for bucket 21
      Entering block accumulator loop for bucket 21:
      bucket 21: 10%
      bucket 20: 20%
      bucket 21: 20%
      bucket 20: 30%
    Getting block 22 of 28
      Reserving size (122828718) for bucket 22
      Calculating Z arrays for bucket 22
      Entering block accumulator loop for bucket 22:
      bucket 21: 30%
      bucket 22: 10%
      bucket 20: 40%
      bucket 21: 40%
      bucket 22: 20%
      bucket 20: 50%
      bucket 21: 50%
      bucket 22: 30%
    Getting block 23 of 28
      Reserving size (122828718) for bucket 23
      Calculating Z arrays for bucket 23
      Entering block accumulator loop for bucket 23:
      bucket 20: 60%
      bucket 21: 60%
      bucket 22: 40%
      bucket 23: 10%
      bucket 21: 70%
      bucket 20: 70%
      bucket 22: 50%
      bucket 23: 20%
      bucket 21: 80%
      bucket 22: 60%
      bucket 20: 80%
      bucket 23: 30%
      bucket 22: 70%
      bucket 21: 90%
      bucket 20: 90%
      bucket 23: 40%
      bucket 22: 80%
      bucket 21: 100%
      Sorting block of length 107637115 for bucket 21
      (Using difference cover)
      bucket 20: 100%
      Sorting block of length 121714038 for bucket 20
      (Using difference cover)
      bucket 23: 50%
      bucket 22: 90%
      bucket 23: 60%
      bucket 22: 100%
      Sorting block of length 45326504 for bucket 22
      (Using difference cover)
      bucket 23: 70%
      bucket 23: 80%
      bucket 23: 90%
      bucket 23: 100%
      Sorting block of length 120387862 for bucket 23
      (Using difference cover)
      Sorting block time: 00:01:04
    Returning block of 45326505 for bucket 22
    Getting block 24 of 28
      Reserving size (122828718) for bucket 24
      Calculating Z arrays for bucket 24
      Entering block accumulator loop for bucket 24:
      bucket 24: 10%
      bucket 24: 20%
      bucket 24: 30%
      bucket 24: 40%
      bucket 24: 50%
      bucket 24: 60%
      bucket 24: 70%
      bucket 24: 80%
      bucket 24: 90%
      bucket 24: 100%
      Sorting block of length 98070715 for bucket 24
      (Using difference cover)
      Sorting block time: 00:02:41
    Returning block of 107637116 for bucket 21
      Sorting block time: 00:03:06
    Returning block of 121714039 for bucket 20
    Getting block 25 of 28
      Reserving size (122828718) for bucket 25
      Calculating Z arrays for bucket 25
      Entering block accumulator loop for bucket 25:
      bucket 25: 10%
      bucket 25: 20%
      bucket 25: 30%
      Sorting block time: 00:03:00
    Returning block of 120387863 for bucket 23
      bucket 25: 40%
      bucket 25: 50%
      bucket 25: 60%
      bucket 25: 70%
    Getting block 26 of 28
      Reserving size (122828718) for bucket 26
      Calculating Z arrays for bucket 26
      Entering block accumulator loop for bucket 26:
      bucket 25: 80%
      bucket 26: 10%
      bucket 25: 90%
      bucket 26: 20%
      bucket 25: 100%
      Sorting block of length 91071917 for bucket 25
      (Using difference cover)
      bucket 26: 30%
    Getting block 27 of 28
      Reserving size (122828718) for bucket 27
      Calculating Z arrays for bucket 27
      Entering block accumulator loop for bucket 27:
      bucket 26: 40%
      bucket 27: 10%
      bucket 26: 50%
      bucket 27: 20%
      bucket 26: 60%
      bucket 27: 30%
      bucket 26: 70%
      bucket 27: 40%
      bucket 26: 80%
      bucket 27: 50%
      bucket 26: 90%
      bucket 27: 60%
      bucket 26: 100%
      Sorting block of length 54558042 for bucket 26
      (Using difference cover)
      bucket 27: 70%
      bucket 27: 80%
      bucket 27: 90%
      bucket 27: 100%
      Sorting block of length 122017027 for bucket 27
      (Using difference cover)
      Sorting block time: 00:02:34
    Returning block of 98070716 for bucket 24
    Getting block 28 of 28
      Reserving size (122828718) for bucket 28
      Calculating Z arrays for bucket 28
      Entering block accumulator loop for bucket 28:
      bucket 28: 10%
      bucket 28: 20%
      bucket 28: 30%
      bucket 28: 40%
      bucket 28: 50%
      bucket 28: 60%
      bucket 28: 70%
      bucket 28: 80%
      bucket 28: 90%
      bucket 28: 100%
      Sorting block of length 96934359 for bucket 28
      (Using difference cover)
      Sorting block time: 00:01:36
    Returning block of 54558043 for bucket 26
      Sorting block time: 00:02:34
    Returning block of 91071918 for bucket 25
      Sorting block time: 00:03:26
    Returning block of 122017028 for bucket 27
      Sorting block time: 00:02:36
    Returning block of 96934360 for bucket 28
    Exited Ebwt loop
    fchr[A]: 0
    fchr[C]: 763450008
    fchr[G]: 1309739408
    fchr[T]: 1856187599
    fchr[$]: 2620345972
    Exiting Ebwt::buildToDisk()
    Returning from initFromVector
    Wrote 877656917 bytes to primary EBWT file: /input_data/mm9.rev.1.bt2
    Wrote 655086500 bytes to secondary EBWT file: /input_data/mm9.rev.2.bt2
    Re-opening _in1 and _in2 as input streams
    Returning from Ebwt constructor
    Headers:
        len: 2620345972
        bwtLen: 2620345973
        sz: 655086493
        bwtSz: 655086494
        lineRate: 6
        offRate: 4
        offMask: 0xfffffff0
        ftabChars: 10
        eftabLen: 20
        eftabSz: 80
        ftabLen: 1048577
        ftabSz: 4194308
        offsLen: 163771624
        offsSz: 655086496
        lineSz: 64
        sideSz: 64
        sideBwtSz: 48
        sideBwtLen: 192
        numSides: 13647636
        numLines: 13647636
        ebwtTotLen: 873448704
        ebwtTotSz: 873448704
        color: 0
        reverse: 1
    Total time for backward call to driver() for mirror index: 00:32:28



```bash
#BWA
bwa index -a bwtsw $fasta_in
```


```bash
#Samtools
samtools faidx reference.fa 
```


```bash
#GATK
java -jar picard.jar CreateSequenceDictionary \
    REFERENCE=$fasta_in \ 
    OUTPUT=reference.dict 
```


```bash
#STAR -- transcriptome
```


```bash
#vg -- genome
```
