```python
from Bio.Seq import Seq
# must go to terminal launcher to first import Biopython by just typing: 
# pip install Biopython
```


```python
my_seq = Seq("GATCG")
```


```python
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G



```python
# We can also print the length of each sequence
print(len(my_seq))
```

    5



```python
# We can also access the different positions of the seq objects 
print(my_seq[0])
```

    G



```python
# Will print the 4th letter in the sequence 
print(my_seq[4])
```

    G



```python
# We only have 4 letters in the sequence so an error will display when we go beyond that 
# ex: print(my_seq[5]); following is a corrected seq 
print(my_seq[2])
```

    T



```python
# Dot count function 
Seq("AAAA").count("AA")
```




    2




```python
# personal random sequence 
my_seq = Seq("GATTGCTTCCCGGATTGCGTACGTACCGAATCCG")
```


```python
# length of seq 
len(my_seq)
```




    34




```python
# count specific letters in a seq 
my_seq.count("G")
```




    9




```python
# Percentage of seq that is G and C 
100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
```




    55.88235294117647




```python
# import gc function to hopefully do that of what is above 
from Bio.SeqUtils import gc_fraction
```


```python
my_seq = Seq("GATTGCTTCCCGGATTGCGTACGTACCGAATCCG")
```


```python
# should get the same value as the percentage above 
gc_fraction(my_seq)
```




    0.5588235294117647




```python
# cut a snip of the seq out and display it 
my_seq[4:12]
```




    Seq('GCTTCCCG')




```python
# slice it with a start, stop, and stride 
# start at place 0, skip 2, then print then repeat 
my_seq[0::3]
```




    Seq('GTTCGTGCAGTG')




```python
my_seq[1::3]
```




    Seq('AGTCAGTGCAC')




```python
my_seq[2:3]
```




    Seq('T')




```python
# we can use negatives to start from the other end 
# this prints the entire thing, but backwards 
my_seq[::-1]
```




    Seq('GCCTAAGCCATGCATGCGTTAGGCCCTTCGTTAG')




```python
# how to turn a seq object back into a string 
str(my_seq)
```




    'GATTGCTTCCCGGATTGCGTACGTACCGAATCCG'




```python
# have a verb placeholder for string formatting (labelling)
fasta_format_string = ">Name\n%s\n" % my_seq
```


```python
print(fasta_format_string)
```

    >Name
    GATTGCTTCCCGGATTGCGTACGTACCGAATCCG
    



```python
# adding two scritps together 
seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
```


```python
seq1 + seq2
```




    Seq('ACGTAACCGG')




```python
seq2 + seq1
```




    Seq('AACCGGACGT')




```python
 #maniulate strings 
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTGCA")]
```


```python
# N is used towards the ends of th strands normally because it gets too small to read 
spacer = Seq("N" * 10)
```


```python
# take the spacer that we made, join with contigs 
# joins in between the seqs
spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTGCA')




```python
dna_seq = Seq("acgtACGT")
```


```python
dna_seq
```




    Seq('acgtACGT')




```python
# makes everything in a sequence become uppercase 
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
# makes everything in a sequence become lowercase 
dna_seq.lower()
```




    Seq('acgtacgt')




```python
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
# search for a certain tyep of codon 
"gtac" in dna_seq
```




    False




```python
# is case sensitive 
"GTAC" in dna_seq
```




    False




```python
dna_seq
```




    Seq('acgtACGT')




```python
# Must save new line as previous name if changing a characteristic of the seq 
dna_seq = dna_seq.upper()
```


```python
"GTAC" in dna_seq
```




    True




```python
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
# creates the complement
my_seq.complement()
```




    Seq('CTAGCTACCCGGATATATCCTAGCTTTTAGCG')




```python
# creates the reverse complement 
my_seq.reverse_complement()
```




    Seq('GCGATTTTCGATCCTATATAGGCCCATCGATC')




```python
# can get the complement for protein seq, is biologically meaningless in reverse 
# gets these from ambiguity code for nucleotides 

protein_seq = Seq("EVRNAK")
protein_seq.complement()

# runs through the seq and realized it could be any of these three if it shows a this one letter 
# V means A C or G 
```




    Seq('EBYNTM')




```python
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
```


```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
template_dna = coding_dna.reverse_complement()
```


```python
# normally how transcription works 
template_dna
```




    Seq('CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT')




```python
# translate into a protein 
messenger_rna = coding_dna.transcribe()
```


```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
template_dna.reverse_complement().transcribe()
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
# turn coding_dna into RNA with transcribe function, back into original DNA with back_transcribe function (reverse transcription)
messenger_rna.back_transcribe()
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
# translate = RNA into protein 
# * indicate a stop codon 
messenger_rna.translate()
```




    Seq('MAIVMGR*KGAR*')




```python
# how to do translation of mitochondrial genome 
coding_dna.translate(table="Vertebrate Mitochondrial")
```




    Seq('MAIVMGRWKGAR*')




```python
# every table is designated a number if you know it, this is a simpilier way to do it 
coding_dna.translate(table = 2)
```




    Seq('MAIVMGRWKGAR*')




```python
# if you just want it to stop at the stop codon 
coding_dna.translate(to_stop = True)
```




    Seq('MAIVMGR')




```python
# you would expect that we'd get the same thing as table 2 above because the sto pcodon is at the end 
coding_dna.translate(table = 2, to_stop = True)
```




    Seq('MAIVMGRWKGAR')




```python
# can create/chnage stop codon to what we want 
coding_dna.translate(table = 2, stop_symbol = "!")
```




    Seq('MAIVMGRWKGAR!')




```python
# what if our sequence uses a non-standard start codon 
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGATAATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACATTATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")
```


```python
gene.translate(table = "Bacterial")
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HR*')




```python
# published all the way to the stop codon
gene.translate(table = "Bacterial", to_stop = True)
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
# telling Biopython that this gene is complete and should start codon with Mythyenine not V, changes it 
gene.translate(table = "Bacterial", cds = True)
```




    Seq('MKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
# codon tables 
from Bio.Data import CodonTable
```


```python
# naming the table for table 1
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```


```python
# nameing the table for table 2 
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
# visualize codon table 
print(standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
# can use tables to remind us of stop and start codons 
mito_table.stop_codons
```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
mito_table.start_codons
```




    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']




```python
# Sequence comparison 
seq = Seq("ACGT")
```


```python
# is fully true 
"ACGT" == seq1
```




    True




```python
# is fully true 
seq1 == "ACGT"
```




    True




```python
# code where length is knoen but not letters 
unknown_seq = Seq(None, 10)
```


```python
unknown_seq
```




    Seq(None, length=10)




```python
# length of unknown variable 
len(unknown_seq)
```




    10




```python
# pulled human gene from file 
seq = Seq({117512683: "TTGAAAACCTGAATGTGCGAGTCAGTCAAGGATAGT"}, length = 159345973)
```


```python
# we don't haveit defined since it's wayyyy upstream
seq[1000:1020]
```




    Seq(None, length=20)




```python
# will display because we have it defined 
seq[117512690:117512700]
```




    Seq('CCTGAATGTG')




```python
# gives us from that point until the end 
seq[117512670:]
```




    Seq({13: 'TTGAAAACCTGAATGTGCGAGTCAGTCAAGGATAGT'}, length=41833303)




```python
seq = Seq("ACGT")
```


```python
undefined_seq = Seq(None, length = 10)
```


```python
# adding seqs together 
seq + undefined_seq +  seq
```




    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)




```python
# we have to make sure we don't change the seq data 
# means it is unmutable 
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
```


```python
# this is on purpose so that throughout the code we aren't able ot accidentally change a codon 
# my_seq[5] = "G"
```


```python
# can create mutable sequence  
from Bio.Seq import MutableSeq
```


```python
mutable_seq = MutableSeq(my_seq)
```


```python
mutable_seq
```




    MutableSeq('GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
mutable_seq[5] = "C"
```


```python
# see that the T in place  5 was changed to a C 
mutable_seq
```




    MutableSeq('GCCATCGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
# can remove the first one of a certain AA from the seq
mutable_seq.remove("T")
```


```python
mutable_seq
```




    MutableSeq('GCCACGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
# can get the reverse of the seq
mutable_seq.reverse()
```


```python
mutable_seq
```




    MutableSeq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
# takes it back to not being able to be changed/ muted 
new_seq = Seq(mutable_seq)
```


```python
new_seq
```




    Seq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
# we imported these for strings in case you don't want to install Biopython built in functions 
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
```


```python
my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG"
```


```python
reverse_complement(my_string)
```




    'CTAACCAGCAGCACGACCACCCTTCCAACGACCCATAACAGC'




```python
transcribe(my_string)
```




    'GCUGUUAUGGGUCGUUGGAAGGGUGGUCGUGCUGCUGGUUAG'




```python
back_transcribe(my_string)
```




    'GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG'




```python
translate(my_string)
```




    'AVMGRWKGGRAAG*'




```python

```
