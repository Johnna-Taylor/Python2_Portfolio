# Python2_Portfolio
This is the portfolio of python that I learned during BISC450C-V84

## Sequence Alignment: Pairwise Alignment 
## Parts 4 & 5

In this class section we learned how to target a codon, target a query, do reverse complements, and search for codon sequences

```python
# before running, go to Terminal and type: pip install Biopython 
from Bio import Align
```


```python
aligner = Align.PairwiseAligner()
```


```python
aligner = Align.PairwiseAligner(match_score = 1.0)
```


```python
target = "GAACT"
```


```python
query = "GAT"
```


```python
score = aligner.score(target,query)
```


```python
score
```




    3.0




```python
alignments = aligner.align(target, query)
```


```python
# we can go in order or skip 
for alignment in alignments:
    print(alignment)
```

    target            0 GAACT 5
                      0 ||--| 5
    query             0 GA--T 3
    
    target            0 GAACT 5
                      0 |-|-| 5
    query             0 G-A-T 3
    



```python
aligner.mode = "local"
```


```python
target = "AGAACTC"
```


```python
query = "GAACT"
```


```python
score = aligner.score(target, query)
```


```python
score
```




    5.0




```python
alignments = aligner.align(target,query)
```


```python
# separates the first and last letter, and isolates the target that we havein the middle shown by lines 
for alignment in alignments:
    print(alignment)
```

    target            1 GAACT 6
                      0 ||||| 5
    query             0 GAACT 5
    



```python
#indels or transposable elements can be found easier with this 
print(aligner)
```

    Pairwise sequence aligner with parameters
      wildcard: None
      match_score: 1.000000
      mismatch_score: 0.000000
      target_internal_open_gap_score: 0.000000
      target_internal_extend_gap_score: 0.000000
      target_left_open_gap_score: 0.000000
      target_left_extend_gap_score: 0.000000
      target_right_open_gap_score: 0.000000
      target_right_extend_gap_score: 0.000000
      query_internal_open_gap_score: 0.000000
      query_internal_extend_gap_score: 0.000000
      query_left_open_gap_score: 0.000000
      query_left_extend_gap_score: 0.000000
      query_right_open_gap_score: 0.000000
      query_right_extend_gap_score: 0.000000
      mode: local
    



```python
aligner.algorithm
```




    'Smith-Waterman'




```python
# significance number 
aligner.epsilon
```




    1e-06




```python
aligner = Align.PairwiseAligner()
```


```python
target = "GAACT"
```


```python
query = "GAT"
```


```python
alignments = aligner.align(target, query)
```


```python
alignment = alignments[0]
```


```python
# like Seq obects, we're looknig at alignment objects 
# beyond our actual Seq
alignment
```




    <Alignment object (2 rows x 5 columns) at 0x7f206960cb90>




```python
alignment.score
```




    3.0




```python
alignment.target
```




    'GAACT'




```python
alignment.query
```




    'GAT'




```python
print(alignment)
```

    target            0 GAACT 5
                      0 ||--| 5
    query             0 GA--T 3
    



```python
alignment.coordinates
```




    array([[0, 2, 4, 5],
           [0, 2, 2, 3]])




```python
# gets the length 
len(alignment)
```




    2




```python
# shape of alignment 
alignment.shape
```




    (2, 5)




```python
aligner.mode = "local"
```


```python
local_alignments = aligner.align("TGAACT", "GAC")
```


```python
local_alignment = local_alignments[0]
```


```python
# it drops off the letters on the end 
print(local_alignment)
```

    target            1 GAAC 5
                      0 ||-| 4
    query             0 GA-C 3
    



```python
# rows by columns
local_alignment.shape
```




    (2, 4)




```python
aligner.mode = "global"
```


```python
aligner = Align.PairwiseAligner(match = 1.0, mismatch_score = -10)
```


```python
print(aligner)
```

    Pairwise sequence aligner with parameters
      wildcard: None
      match_score: 1.000000
      mismatch_score: -10.000000
      target_internal_open_gap_score: 0.000000
      target_internal_extend_gap_score: 0.000000
      target_left_open_gap_score: 0.000000
      target_left_extend_gap_score: 0.000000
      target_right_open_gap_score: 0.000000
      target_right_extend_gap_score: 0.000000
      query_internal_open_gap_score: 0.000000
      query_internal_extend_gap_score: 0.000000
      query_left_open_gap_score: 0.000000
      query_left_extend_gap_score: 0.000000
      query_right_open_gap_score: 0.000000
      query_right_extend_gap_score: 0.000000
      mode: global
    



```python
alignments = aligner.align("AAACAAA", "AAAGAAA")
```


```python
len(alignment)
```




    2




```python
# score is 0 would rather put in a mismatch than a zero 
print(alignments[0])
```

    target            0 AAAC-AAA 7
                      0 |||--||| 8
    query             0 AAA-GAAA 7
    



```python
# reverses G and C 
print(alignments[1])
```

    target            0 AAA-CAAA 7
                      0 |||--||| 8
    query             0 AAAG-AAA 7
    



```python
print(local_alignment)
```

    target            1 GAAC 5
                      0 ||-| 4
    query             0 GA-C 3
    



```python
local_alignment.sort()
```


```python
print(local_alignment)
```

    target            0 GA-C 3
                      0 ||-| 4
    query             1 GAAC 5
    



```python
from Bio.Seq import reverse_complement
```


```python
target = "AAACCC"
```


```python
query = "AACC"
```


```python
aligner = Align.PairwiseAligner(mismatch_score = -1, internal_gap_score = -1)
```


```python
# align to positive strand 
aligner.score(target, query)
```




    4.0




```python
# trying to align to reverse complement - negative strand 
aligner.score(target, reverse_complement(query))
```




    0.0




```python
# aligning twwo negatives to make a positive 
aligner.score(target, reverse_complement(query), strand = "-")
```




    4.0




```python
aligner.score(target, query, strand = "-")
```




    0.0




```python
alignments = aligner.align(target, query)
```


```python
len(alignments)
```




    1




```python
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             0 -AACC- 4
    



```python
# can save to a certian type of file 
print(alignments[0].format("bed") )
```

    target	1	5	query	4.0	+	1	5	0	1	4,	0,
    



```python
alignments = aligner.align(target, reverse_complement(query), strand = "-")
```


```python
print(alignments[0].format("bed"))
```

    target	1	5	query	4.0	-	1	5	0	1	4,	0,
    



```python
alignments = aligner.align(target, query, strand = "-")
```


```python
len(alignments)
```




    2




```python
print(alignments[0])
```

    target            0 AAACCC----  6
                      0 ---------- 10
    query             4 ------GGTT  0
    



```python
# there is a super negative score, so there is zero alignment 
print(alignments[1])
```

    target            0 ----AAACCC  6
                      0 ---------- 10
    query             4 GGTT------  0
    



```python
aligner.left_gap_score = -0.5
```


```python
aligner.right_gap_score = -0.2
```


```python
aligner.score(target,query)
```




    3.3




```python
alignments = aligner.align(target, query)
```


```python
len(alignments)
```




    1




```python
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             0 -AACC- 4
    



```python
alignments = aligner.align(target, reverse_complement(query), strand = "-")
```


```python
print(alignments)
```

    <Bio.Align.PairwiseAlignments object at 0x7f20695de7d0>



```python
aligner.score(target, reverse_complement(query), strand = "-")
```




    3.3




```python
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             4 -AACC- 0
    



```python
# shows that mapping/printing our DNa gives us a worse alignment than mapping to our standard string 
aligner.score(target, reverse_complement(query), strand = "+")
```




    -2.0




```python
# will give you the correct alignment 
aligner.score(target, query, strand = "+")
```




    3.3


