

```python
input_dir = "/coop_SE/primed_testing/"
output_dir = "/coop_SE/primed_testing/"
input_fq_1 = "ATAC_h7_Smtmrs_1_1M.fq.gz"
input_fq_2 = "ATAC_h7_Smtmrs_2_1M.fq.gz"
```


```python
!zcat $input_dir$input_fq_2 | head
```

    @SRR3689991.1.2 J00118:114:H3WLKBBXX:4:1101:2260:1226 length=151
    NCACGGGCTCGGGCNAGGCTTGCTGANANTCTCAGGGCACGGGGGATCCAGGAGTCCAAGGTGGACCCAGGGGGCCTGGCCGCCTGCACAGCCCCTTNTTGTCNTNGATNGGAAGGAGCCCNTCAGACCCCTGTCCCNTAGNCACCTCCGA
    +SRR3689991.1.2 J00118:114:H3WLKBBXX:4:1101:2260:1226 length=151
    #AAAF--A<FFJ<J#JJFAFF7F-7<#J#-FF7J--FF-<A-FAJJFFFJ---7-7J<77AAJ<AJ-AAF<7A-77<JFFJJAJA-7FAJFJF77-7#A<<-F#J#A--#7---77AF-7-#-7F---FFF)))7)F#---#7----7)-)
    @SRR3689991.2.2 J00118:114:H3WLKBBXX:4:1101:3315:1226 length=151
    NGTTGATTCGGGNGNATCCTATTGGTNCNGGGGCTTTGTATGATTATGGGCGGTGATTAGTAGTAGTTACCTGTCTCTTATACACAATCGACGCTTCNGACGANTNTGGNTCTCGGTTCTCNTCGTATCATCGATAANAGANGCCTCATTC
    +SRR3689991.2.2 J00118:114:H3WLKBBXX:4:1101:3315:1226 length=151
    #-AA-<FJJAFA#F#FJ7<FJJF-FJ#-#7-AF<-7F--FF7JAJFJFJF-A---7A<F--77AF-7AF---<-7A<FF<JFFA--<--77AF-7--#7AF77#-#7--#7A77-----7-#--<F7F7-7----7-#--7#--)-7)--7
    @SRR3689991.3.2 J00118:114:H3WLKBBXX:4:1101:3356:1226 length=151
    NGTCTGTGGCGANCNCGATTCTCCCTNTNGGGTGGCTACAGGCTAGAAACTGTCTCTTATACACATCTGACGCTGCCGACGAGTGTAGATCTCGGTGNTCGCCNTNTCANTAGAAAAAAAANCTCGTACCCGCGCCCNCTCNGGACCCCTC
    
    gzip: stdout: Broken pipe



```python
!cutadapt \
            -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
            -o $input_dir"trimmed_"$input_fq_1 -p $input_dir"trimmed_"$input_fq_2 \
            $input_dir$input_fq_1 $input_dir$input_fq_2
```

    This is cutadapt 1.13 with Python 2.7.13
    Command line parameters: -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o /coop_SE/primed_testing/trimmed_ATAC_h7_Smtmrs_1_1M.fq.gz -p /coop_SE/primed_testing/trimmed_ATAC_h7_Smtmrs_2_1M.fq.gz /coop_SE/primed_testing/ATAC_h7_Smtmrs_1_1M.fq.gz /coop_SE/primed_testing/ATAC_h7_Smtmrs_2_1M.fq.gz
    Trimming 2 adapters with at most 10.0% errors in paired-end mode ...
    Finished in 45.29 s (45 us/read; 1.32 M reads/minute).
    
    === Summary ===
    
    Total read pairs processed:          1,000,000
      Read 1 with adapter:                  23,518 (2.4%)
      Read 2 with adapter:                  30,681 (3.1%)
    Pairs written (passing filters):     1,000,000 (100.0%)
    
    Total basepairs processed:   302,000,000 bp
      Read 1:   151,000,000 bp
      Read 2:   151,000,000 bp
    Total written (filtered):    301,805,485 bp (99.9%)
      Read 1:   150,916,355 bp
      Read 2:   150,889,130 bp
    
    === First read: Adapter 1 ===
    
    Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 23518 times.
    
    No. of allowed errors:
    0-9 bp: 0; 10-13 bp: 1
    
    Bases preceding removed adapters:
      A: 34.7%
      C: 21.7%
      G: 30.0%
      T: 13.6%
      none/other: 0.0%
    
    Overview of removed sequences
    length	count	expect	max.err	error counts
    3	20214	15625.0	0	20214
    4	2345	3906.2	0	2345
    5	581	976.6	0	581
    6	73	244.1	0	73
    7	16	61.0	0	16
    8	6	15.3	0	6
    9	18	3.8	0	11 7
    10	31	1.0	1	7 24
    11	9	0.2	1	4 5
    12	11	0.1	1	6 5
    13	3	0.0	1	3
    14	9	0.0	1	5 4
    15	11	0.0	1	8 3
    16	6	0.0	1	4 2
    17	6	0.0	1	6
    18	6	0.0	1	2 4
    19	8	0.0	1	7 1
    20	6	0.0	1	3 3
    21	5	0.0	1	5
    22	5	0.0	1	0 5
    23	2	0.0	1	1 1
    24	4	0.0	1	3 1
    25	4	0.0	1	4
    26	5	0.0	1	3 2
    27	7	0.0	1	5 2
    28	10	0.0	1	9 1
    29	4	0.0	1	4
    30	1	0.0	1	0 1
    31	3	0.0	1	2 1
    32	3	0.0	1	2 1
    33	2	0.0	1	2
    34	5	0.0	1	4 1
    35	6	0.0	1	6
    36	1	0.0	1	0 1
    37	2	0.0	1	2
    38	4	0.0	1	4
    39	3	0.0	1	3
    40	2	0.0	1	2
    41	4	0.0	1	4
    42	4	0.0	1	4
    43	4	0.0	1	2 2
    44	2	0.0	1	2
    45	3	0.0	1	3
    46	3	0.0	1	3
    49	1	0.0	1	1
    50	1	0.0	1	1
    52	1	0.0	1	1
    53	4	0.0	1	4
    54	2	0.0	1	1 1
    55	1	0.0	1	1
    56	1	0.0	1	0 1
    57	2	0.0	1	2
    59	3	0.0	1	3
    60	2	0.0	1	2
    62	2	0.0	1	0 2
    64	1	0.0	1	1
    65	1	0.0	1	1
    67	1	0.0	1	1
    69	1	0.0	1	1
    70	3	0.0	1	2 1
    71	2	0.0	1	1 1
    76	1	0.0	1	1
    85	1	0.0	1	1
    87	2	0.0	1	0 2
    91	1	0.0	1	1
    93	1	0.0	1	0 1
    96	1	0.0	1	0 1
    98	2	0.0	1	0 2
    99	3	0.0	1	1 2
    100	1	0.0	1	1
    104	2	0.0	1	0 2
    105	1	0.0	1	0 1
    106	1	0.0	1	0 1
    108	1	0.0	1	1
    115	1	0.0	1	1
    120	1	0.0	1	0 1
    122	1	0.0	1	0 1
    124	1	0.0	1	0 1
    128	1	0.0	1	0 1
    130	2	0.0	1	0 2
    138	2	0.0	1	0 2
    142	1	0.0	1	0 1
    144	1	0.0	1	0 1
    146	1	0.0	1	0 1
    151	2	0.0	1	1 1
    
    === Second read: Adapter 2 ===
    
    Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 30681 times.
    
    No. of allowed errors:
    0-9 bp: 0; 10-13 bp: 1
    
    Bases preceding removed adapters:
      A: 34.0%
      C: 18.8%
      G: 10.7%
      T: 36.5%
      none/other: 0.0%
    
    Overview of removed sequences
    length	count	expect	max.err	error counts
    3	22562	15625.0	0	22562
    4	4574	3906.2	0	4574
    5	3190	976.6	0	3190
    6	79	244.1	0	79
    7	30	61.0	0	30
    8	4	15.3	0	4
    9	19	3.8	0	5 14
    10	29	1.0	1	6 23
    11	11	0.2	1	2 9
    12	11	0.1	1	4 7
    13	5	0.0	1	3 2
    14	9	0.0	1	3 6
    15	6	0.0	1	5 1
    16	4	0.0	1	1 3
    17	4	0.0	1	3 1
    18	3	0.0	1	1 2
    19	6	0.0	1	1 5
    20	5	0.0	1	3 2
    21	4	0.0	1	1 3
    22	1	0.0	1	0 1
    23	2	0.0	1	1 1
    24	3	0.0	1	2 1
    25	4	0.0	1	3 1
    26	3	0.0	1	2 1
    27	4	0.0	1	3 1
    28	7	0.0	1	6 1
    29	3	0.0	1	3
    31	2	0.0	1	1 1
    32	4	0.0	1	1 3
    33	2	0.0	1	2
    34	1	0.0	1	1
    35	4	0.0	1	4
    36	1	0.0	1	1
    37	3	0.0	1	2 1
    38	6	0.0	1	3 3
    39	2	0.0	1	2
    40	3	0.0	1	2 1
    41	4	0.0	1	4
    42	3	0.0	1	2 1
    43	3	0.0	1	3
    44	3	0.0	1	2 1
    45	3	0.0	1	3
    46	3	0.0	1	2 1
    49	1	0.0	1	1
    50	1	0.0	1	1
    52	1	0.0	1	1
    53	4	0.0	1	3 1
    54	2	0.0	1	1 1
    55	1	0.0	1	1
    56	1	0.0	1	1
    57	2	0.0	1	2
    58	1	0.0	1	0 1
    59	3	0.0	1	3
    60	3	0.0	1	1 2
    64	1	0.0	1	1
    65	1	0.0	1	1
    66	1	0.0	1	1
    67	1	0.0	1	1
    69	2	0.0	1	1 1
    70	4	0.0	1	3 1
    71	1	0.0	1	1
    76	1	0.0	1	1
    77	2	0.0	1	0 2
    85	1	0.0	1	1
    91	1	0.0	1	1
    93	2	0.0	1	1 1
    99	1	0.0	1	1
    100	1	0.0	1	1
    106	1	0.0	1	0 1
    108	1	0.0	1	1
    111	1	0.0	1	0 1
    115	1	0.0	1	1
    128	1	0.0	1	0 1
    130	1	0.0	1	0 1
    131	1	0.0	1	0 1
    139	2	0.0	1	0 2
    146	1	0.0	1	0 1
    150	1	0.0	1	0 1
    151	1	0.0	1	1
    



```python
!/root/sickle/sickle pe -f $input_dir"trimmed_"$input_fq_1 -r $input_dir"trimmed_"$input_fq_2 \
    -o $input_dir"filt_"$input_fq_1 -p $input_dir"filt_"$input_fq_2 -s $input_dir"single_"$input_fq_1\
    -t sanger -q 20 -l 20 -g
```

    
    PE forward file: /coop_SE/primed_testing/trimmed_ATAC_h7_Smtmrs_1_1M.fq.gz
    PE reverse file: /coop_SE/primed_testing/trimmed_ATAC_h7_Smtmrs_2_1M.fq.gz
    
    Total input FastQ records: 2000000 (1000000 pairs)
    
    FastQ paired records kept: 1959084 (979542 pairs)
    FastQ single records kept: 20400 (from PE1: 20284, from PE2: 116)
    FastQ paired records discarded: 116 (58 pairs)
    FastQ single records discarded: 20400 (from PE1: 116, from PE2: 20284)
    



```python
!zcat $input_dir"filt_"$input_fq_1 | tail
```

    +
    A<AF<FA-FAJJJF7FFJJJJJJJJ7AJ77-AFJA-A-<<FA--<FFJJJ<JJJF--<AJJF<AAJA<-<-AAJF77<7FJFJJFJ<JF7JF7FJ--77FF--<<7AF<<A-FAAA-7FJJJJJ777<J<A<FFA7
    @SRR3689991.999999.1 J00118:114:H3WLKBBXX:4:1102:13372:3178 length=151
    TGTTCAAACTGTCATTTTATTTTTACGTTGTTAGATATGGGGAGTAGTGTGATTGCTGTCTCTTATACACATCTCCGAGCCCACGAGACCTCTCTACATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAA
    +
    AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJFFJJJJJJJJJJJJJJJJJJJFJF-FJJJJJJFJJAAFFJ<AFFFJAAFJ<AFJJF7
    @SRR3689991.1000000.1 J00118:114:H3WLKBBXX:4:1102:13433:3178 length=151
    GACTAATCTTCAACTCCTACATACTTCCCCCATTATTCCTAGAACCTGTCTCTTATACACATCTCCGAGCCCACGAGACCTCTCTACATCTCGTATGCCGTCTTCTGCTTGAAAAAAAA
    +
    AAFFFFJJJJJJJJJJJJJJJJJJJJJJJJJJA<FFFJJJJFJFJJJJ<A7FJFJJJFJJJJJJJJJF<JJJJJFJ<JFAFFJJJJJJ-7AFJJJJJAJJJF-AFFJFFFJJ<<JJJJJ



```python

```
