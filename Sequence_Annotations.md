# Python2_Portfolio
This is the portfolio of python that I learned during BISC450C-V84

## Sequence Annotations
### Parts 1, 2, 3, & 4

In this class we discussed how to do Sequence Annotations by creating a seq record, annotating it, and displaying some data from it

```python
# make sure to go to terminal and do: pip install Biopython first to run this 
from Bio.SeqRecord import SeqRecord
```


```python
help(SeqRecord)
```

    Help on class SeqRecord in module Bio.SeqRecord:
    
    class SeqRecord(builtins.object)
     |  SeqRecord(seq, id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=None, features=None, annotations=None, letter_annotations=None)
     |  
     |  A SeqRecord object holds a sequence and information about it.
     |  
     |  Main attributes:
     |   - id          - Identifier such as a locus tag (string)
     |   - seq         - The sequence itself (Seq object or similar)
     |  
     |  Additional attributes:
     |   - name        - Sequence name, e.g. gene name (string)
     |   - description - Additional text (string)
     |   - dbxrefs     - List of database cross references (list of strings)
     |   - features    - Any (sub)features defined (list of SeqFeature objects)
     |   - annotations - Further information about the whole sequence (dictionary).
     |     Most entries are strings, or lists of strings.
     |   - letter_annotations - Per letter/symbol annotation (restricted
     |     dictionary). This holds Python sequences (lists, strings
     |     or tuples) whose length matches that of the sequence.
     |     A typical use would be to hold a list of integers
     |     representing sequencing quality scores, or a string
     |     representing the secondary structure.
     |  
     |  You will typically use Bio.SeqIO to read in sequences from files as
     |  SeqRecord objects.  However, you may want to create your own SeqRecord
     |  objects directly (see the __init__ method for further details):
     |  
     |  >>> from Bio.Seq import Seq
     |  >>> from Bio.SeqRecord import SeqRecord
     |  >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
     |  ...                    id="YP_025292.1", name="HokC",
     |  ...                    description="toxic membrane protein")
     |  >>> print(record)
     |  ID: YP_025292.1
     |  Name: HokC
     |  Description: toxic membrane protein
     |  Number of features: 0
     |  Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF')
     |  
     |  If you want to save SeqRecord objects to a sequence file, use Bio.SeqIO
     |  for this.  For the special case where you want the SeqRecord turned into
     |  a string in a particular file format there is a format method which uses
     |  Bio.SeqIO internally:
     |  
     |  >>> print(record.format("fasta"))
     |  >YP_025292.1 toxic membrane protein
     |  MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
     |  <BLANKLINE>
     |  
     |  You can also do things like slicing a SeqRecord, checking its length, etc
     |  
     |  >>> len(record)
     |  44
     |  >>> edited = record[:10] + record[11:]
     |  >>> print(edited.seq)
     |  MKQHKAMIVAIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
     |  >>> print(record.seq)
     |  MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
     |  
     |  Methods defined here:
     |  
     |  __add__(self, other)
     |      Add another sequence or string to this sequence.
     |      
     |      The other sequence can be a SeqRecord object, a Seq object (or
     |      similar, e.g. a MutableSeq) or a plain Python string. If you add
     |      a plain string or a Seq (like) object, the new SeqRecord will simply
     |      have this appended to the existing data. However, any per letter
     |      annotation will be lost:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
     |      >>> print("%s %s" % (record.id, record.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(record.letter_annotations))
     |      ['solexa_quality']
     |      
     |      >>> new = record + "ACT"
     |      >>> print("%s %s" % (new.id, new.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNNACT
     |      >>> print(list(new.letter_annotations))
     |      []
     |      
     |      The new record will attempt to combine the annotation, but for any
     |      ambiguities (e.g. different names) it defaults to omitting that
     |      annotation.
     |      
     |      >>> from Bio import SeqIO
     |      >>> with open("GenBank/pBAD30.gb") as handle:
     |      ...     plasmid = SeqIO.read(handle, "gb")
     |      >>> print("%s %i" % (plasmid.id, len(plasmid)))
     |      pBAD30 4923
     |      
     |      Now let's cut the plasmid into two pieces, and join them back up the
     |      other way round (i.e. shift the starting point on this plasmid, have
     |      a look at the annotated features in the original file to see why this
     |      particular split point might make sense):
     |      
     |      >>> left = plasmid[:3765]
     |      >>> right = plasmid[3765:]
     |      >>> new = right + left
     |      >>> print("%s %i" % (new.id, len(new)))
     |      pBAD30 4923
     |      >>> str(new.seq) == str(right.seq + left.seq)
     |      True
     |      >>> len(new.features) == len(left.features) + len(right.features)
     |      True
     |      
     |      When we add the left and right SeqRecord objects, their annotation
     |      is all consistent, so it is all conserved in the new SeqRecord:
     |      
     |      >>> new.id == left.id == right.id == plasmid.id
     |      True
     |      >>> new.name == left.name == right.name == plasmid.name
     |      True
     |      >>> new.description == plasmid.description
     |      True
     |      >>> new.annotations == left.annotations == right.annotations
     |      True
     |      >>> new.letter_annotations == plasmid.letter_annotations
     |      True
     |      >>> new.dbxrefs == left.dbxrefs == right.dbxrefs
     |      True
     |      
     |      However, we should point out that when we sliced the SeqRecord,
     |      any annotations dictionary or dbxrefs list entries were lost.
     |      You can explicitly copy them like this:
     |      
     |      >>> new.annotations = plasmid.annotations.copy()
     |      >>> new.dbxrefs = plasmid.dbxrefs[:]
     |  
     |  __bool__(self)
     |      Boolean value of an instance of this class (True).
     |      
     |      This behaviour is for backwards compatibility, since until the
     |      __len__ method was added, a SeqRecord always evaluated as True.
     |      
     |      Note that in comparison, a Seq object will evaluate to False if it
     |      has a zero length sequence.
     |      
     |      WARNING: The SeqRecord may in future evaluate to False when its
     |      sequence is of zero length (in order to better match the Seq
     |      object behaviour)!
     |  
     |  __bytes__(self)
     |  
     |  __contains__(self, char)
     |      Implement the 'in' keyword, searches the sequence.
     |      
     |      e.g.
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Fasta/sweetpea.nu", "fasta")
     |      >>> "GAATTC" in record
     |      False
     |      >>> "AAA" in record
     |      True
     |      
     |      This essentially acts as a proxy for using "in" on the sequence:
     |      
     |      >>> "GAATTC" in record.seq
     |      False
     |      >>> "AAA" in record.seq
     |      True
     |      
     |      Note that you can also use Seq objects as the query,
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> Seq("AAA") in record
     |      True
     |      
     |      See also the Seq object's __contains__ method.
     |  
     |  __eq__(self, other)
     |      Define the equal-to operand (not implemented).
     |  
     |  __format__(self, format_spec)
     |      Return the record as a string in the specified file format.
     |      
     |      This method supports the Python format() function and f-strings.
     |      The format_spec should be a lower case string supported by
     |      Bio.SeqIO as a text output file format. Requesting a binary file
     |      format raises a ValueError. e.g.
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
     |      ...                    id="YP_025292.1", name="HokC",
     |      ...                    description="toxic membrane protein")
     |      ...
     |      >>> format(record, "fasta")
     |      '>YP_025292.1 toxic membrane protein\nMKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF\n'
     |      >>> print(f"Here is {record.id} in FASTA format:\n{record:fasta}")
     |      Here is YP_025292.1 in FASTA format:
     |      >YP_025292.1 toxic membrane protein
     |      MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
     |      <BLANKLINE>
     |      
     |      See also the SeqRecord's format() method.
     |  
     |  __ge__(self, other)
     |      Define the greater-than-or-equal-to operand (not implemented).
     |  
     |  __getitem__(self, index)
     |      Return a sub-sequence or an individual letter.
     |      
     |      Slicing, e.g. my_record[5:10], returns a new SeqRecord for
     |      that sub-sequence with some annotation preserved as follows:
     |      
     |      * The name, id and description are kept as-is.
     |      * Any per-letter-annotations are sliced to match the requested
     |        sub-sequence.
     |      * Unless a stride is used, all those features which fall fully
     |        within the subsequence are included (with their locations
     |        adjusted accordingly). If you want to preserve any truncated
     |        features (e.g. GenBank/EMBL source features), you must
     |        explicitly add them to the new SeqRecord yourself.
     |      * With the exception of any molecule type, the annotations
     |        dictionary and the dbxrefs list are not used for the new
     |        SeqRecord, as in general they may not apply to the
     |        subsequence. If you want to preserve them, you must explicitly
     |        copy them to the new SeqRecord yourself.
     |      
     |      Using an integer index, e.g. my_record[5] is shorthand for
     |      extracting that letter from the sequence, my_record.seq[5].
     |      
     |      For example, consider this short protein and its secondary
     |      structure as encoded by the PDB (e.g. H for alpha helices),
     |      plus a simple feature for its histidine self phosphorylation
     |      site:
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> from Bio.SeqFeature import SeqFeature, SimpleLocation
     |      >>> rec = SeqRecord(Seq("MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLAT"
     |      ...                     "EMMSEQDGYLAESINKDIEECNAIIEQFIDYLR"),
     |      ...                 id="1JOY", name="EnvZ",
     |      ...                 description="Homodimeric domain of EnvZ from E. coli")
     |      >>> rec.letter_annotations["secondary_structure"] = "  S  SSSSSSHHHHHTTTHHHHHHHHHHHHHHHHHHHHHHTHHHHHHHHHHHHHHHHHHHHHTT  "
     |      >>> rec.features.append(SeqFeature(SimpleLocation(20, 21),
     |      ...                     type = "Site"))
     |      
     |      Now let's have a quick look at the full record,
     |      
     |      >>> print(rec)
     |      ID: 1JOY
     |      Name: EnvZ
     |      Description: Homodimeric domain of EnvZ from E. coli
     |      Number of features: 1
     |      Per letter annotation for: secondary_structure
     |      Seq('MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLATEMMSEQDGYLAESINKDIEE...YLR')
     |      >>> rec.letter_annotations["secondary_structure"]
     |      '  S  SSSSSSHHHHHTTTHHHHHHHHHHHHHHHHHHHHHHTHHHHHHHHHHHHHHHHHHHHHTT  '
     |      >>> print(rec.features[0].location)
     |      [20:21]
     |      
     |      Now let's take a sub sequence, here chosen as the first (fractured)
     |      alpha helix which includes the histidine phosphorylation site:
     |      
     |      >>> sub = rec[11:41]
     |      >>> print(sub)
     |      ID: 1JOY
     |      Name: EnvZ
     |      Description: Homodimeric domain of EnvZ from E. coli
     |      Number of features: 1
     |      Per letter annotation for: secondary_structure
     |      Seq('RTLLMAGVSHDLRTPLTRIRLATEMMSEQD')
     |      >>> sub.letter_annotations["secondary_structure"]
     |      'HHHHHTTTHHHHHHHHHHHHHHHHHHHHHH'
     |      >>> print(sub.features[0].location)
     |      [9:10]
     |      
     |      You can also of course omit the start or end values, for
     |      example to get the first ten letters only:
     |      
     |      >>> print(rec[:10])
     |      ID: 1JOY
     |      Name: EnvZ
     |      Description: Homodimeric domain of EnvZ from E. coli
     |      Number of features: 0
     |      Per letter annotation for: secondary_structure
     |      Seq('MAAGVKQLAD')
     |      
     |      Or for the last ten letters:
     |      
     |      >>> print(rec[-10:])
     |      ID: 1JOY
     |      Name: EnvZ
     |      Description: Homodimeric domain of EnvZ from E. coli
     |      Number of features: 0
     |      Per letter annotation for: secondary_structure
     |      Seq('IIEQFIDYLR')
     |      
     |      If you omit both, then you get a copy of the original record (although
     |      lacking the annotations and dbxrefs):
     |      
     |      >>> print(rec[:])
     |      ID: 1JOY
     |      Name: EnvZ
     |      Description: Homodimeric domain of EnvZ from E. coli
     |      Number of features: 1
     |      Per letter annotation for: secondary_structure
     |      Seq('MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLATEMMSEQDGYLAESINKDIEE...YLR')
     |      
     |      Finally, indexing with a simple integer is shorthand for pulling out
     |      that letter from the sequence directly:
     |      
     |      >>> rec[5]
     |      'K'
     |      >>> rec.seq[5]
     |      'K'
     |  
     |  __gt__(self, other)
     |      Define the greater-than operand (not implemented).
     |  
     |  __init__(self, seq, id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=None, features=None, annotations=None, letter_annotations=None)
     |      Create a SeqRecord.
     |      
     |      Arguments:
     |       - seq         - Sequence, required (Seq or MutableSeq)
     |       - id          - Sequence identifier, recommended (string)
     |       - name        - Sequence name, optional (string)
     |       - description - Sequence description, optional (string)
     |       - dbxrefs     - Database cross references, optional (list of strings)
     |       - features    - Any (sub)features, optional (list of SeqFeature objects)
     |       - annotations - Dictionary of annotations for the whole sequence
     |       - letter_annotations - Dictionary of per-letter-annotations, values
     |         should be strings, list or tuples of the same length as the full
     |         sequence.
     |      
     |      You will typically use Bio.SeqIO to read in sequences from files as
     |      SeqRecord objects.  However, you may want to create your own SeqRecord
     |      objects directly.
     |      
     |      Note that while an id is optional, we strongly recommend you supply a
     |      unique id string for each record.  This is especially important
     |      if you wish to write your sequences to a file.
     |      
     |      You can create a 'blank' SeqRecord object, and then populate the
     |      attributes later.
     |  
     |  __iter__(self)
     |      Iterate over the letters in the sequence.
     |      
     |      For example, using Bio.SeqIO to read in a protein FASTA file:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Fasta/loveliesbleeding.pro", "fasta")
     |      >>> for amino in record:
     |      ...     print(amino)
     |      ...     if amino == "L": break
     |      X
     |      A
     |      G
     |      L
     |      >>> print(record.seq[3])
     |      L
     |      
     |      This is just a shortcut for iterating over the sequence directly:
     |      
     |      >>> for amino in record.seq:
     |      ...     print(amino)
     |      ...     if amino == "L": break
     |      X
     |      A
     |      G
     |      L
     |      >>> print(record.seq[3])
     |      L
     |      
     |      Note that this does not facilitate iteration together with any
     |      per-letter-annotation.  However, you can achieve that using the
     |      python zip function on the record (or its sequence) and the relevant
     |      per-letter-annotation:
     |      
     |      >>> from Bio import SeqIO
     |      >>> rec = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
     |      >>> print("%s %s" % (rec.id, rec.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(rec.letter_annotations))
     |      ['solexa_quality']
     |      >>> for nuc, qual in zip(rec, rec.letter_annotations["solexa_quality"]):
     |      ...     if qual > 35:
     |      ...         print("%s %i" % (nuc, qual))
     |      A 40
     |      C 39
     |      G 38
     |      T 37
     |      A 36
     |      
     |      You may agree that using zip(rec.seq, ...) is more explicit than using
     |      zip(rec, ...) as shown above.
     |  
     |  __le__(self, other)
     |      Define the less-than-or-equal-to operand (not implemented).
     |  
     |  __len__(self)
     |      Return the length of the sequence.
     |      
     |      For example, using Bio.SeqIO to read in a FASTA nucleotide file:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Fasta/sweetpea.nu", "fasta")
     |      >>> len(record)
     |      309
     |      >>> len(record.seq)
     |      309
     |  
     |  __lt__(self, other)
     |      Define the less-than operand (not implemented).
     |  
     |  __ne__(self, other)
     |      Define the not-equal-to operand (not implemented).
     |  
     |  __radd__(self, other)
     |      Add another sequence or string to this sequence (from the left).
     |      
     |      This method handles adding a Seq object (or similar, e.g. MutableSeq)
     |      or a plain Python string (on the left) to a SeqRecord (on the right).
     |      See the __add__ method for more details, but for example:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
     |      >>> print("%s %s" % (record.id, record.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(record.letter_annotations))
     |      ['solexa_quality']
     |      
     |      >>> new = "ACT" + record
     |      >>> print("%s %s" % (new.id, new.seq))
     |      slxa_0001_1_0001_01 ACTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(new.letter_annotations))
     |      []
     |  
     |  __repr__(self)
     |      Return a concise summary of the record for debugging (string).
     |      
     |      The python built in function repr works by calling the object's __repr__
     |      method.  e.g.
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> rec = SeqRecord(Seq("MASRGVNKVILVGNLGQDPEVRYMPNGGAVANITLATSESWRDKAT"
     |      ...                     "GEMKEQTEWHRVVLFGKLAEVASEYLRKGSQVYIEGQLRTRKWTDQ"
     |      ...                     "SGQDRYTTEVVVNVGGTMQMLGGRQGGGAPAGGNIGGGQPQGGWGQ"
     |      ...                     "PQQPQGGNQFSGGAQSRPQQSAPAAPSNEPPMDFDDDIPF"),
     |      ...                 id="NP_418483.1", name="b4059",
     |      ...                 description="ssDNA-binding protein",
     |      ...                 dbxrefs=["ASAP:13298", "GI:16131885", "GeneID:948570"])
     |      >>> print(repr(rec))
     |      SeqRecord(seq=Seq('MASRGVNKVILVGNLGQDPEVRYMPNGGAVANITLATSESWRDKATGEMKEQTE...IPF'), id='NP_418483.1', name='b4059', description='ssDNA-binding protein', dbxrefs=['ASAP:13298', 'GI:16131885', 'GeneID:948570'])
     |      
     |      At the python prompt you can also use this shorthand:
     |      
     |      >>> rec
     |      SeqRecord(seq=Seq('MASRGVNKVILVGNLGQDPEVRYMPNGGAVANITLATSESWRDKATGEMKEQTE...IPF'), id='NP_418483.1', name='b4059', description='ssDNA-binding protein', dbxrefs=['ASAP:13298', 'GI:16131885', 'GeneID:948570'])
     |      
     |      Note that long sequences are shown truncated. Also note that any
     |      annotations, letter_annotations and features are not shown (as they
     |      would lead to a very long string).
     |  
     |  __str__(self)
     |      Return a human readable summary of the record and its annotation (string).
     |      
     |      The python built in function str works by calling the object's __str__
     |      method.  e.g.
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
     |      ...                    id="YP_025292.1", name="HokC",
     |      ...                    description="toxic membrane protein, small")
     |      >>> print(str(record))
     |      ID: YP_025292.1
     |      Name: HokC
     |      Description: toxic membrane protein, small
     |      Number of features: 0
     |      Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF')
     |      
     |      In this example you don't actually need to call str explicitly, as the
     |      print command does this automatically:
     |      
     |      >>> print(record)
     |      ID: YP_025292.1
     |      Name: HokC
     |      Description: toxic membrane protein, small
     |      Number of features: 0
     |      Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF')
     |      
     |      Note that long sequences are shown truncated.
     |  
     |  count(self, sub, start=None, end=None)
     |      Return the number of non-overlapping occurrences of sub in seq[start:end].
     |      
     |      Optional arguments start and end are interpreted as in slice notation.
     |      This method behaves as the count method of Python strings.
     |  
     |  format(self, format)
     |      Return the record as a string in the specified file format.
     |      
     |      The format should be a lower case string supported as an output
     |      format by Bio.SeqIO, which is used to turn the SeqRecord into a
     |      string.  e.g.
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
     |      ...                    id="YP_025292.1", name="HokC",
     |      ...                    description="toxic membrane protein")
     |      >>> record.format("fasta")
     |      '>YP_025292.1 toxic membrane protein\nMKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF\n'
     |      >>> print(record.format("fasta"))
     |      >YP_025292.1 toxic membrane protein
     |      MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
     |      <BLANKLINE>
     |      
     |      The Python print function automatically appends a new line, meaning
     |      in this example a blank line is shown.  If you look at the string
     |      representation you can see there is a trailing new line (shown as
     |      slash n) which is important when writing to a file or if
     |      concatenating multiple sequence strings together.
     |      
     |      Note that this method will NOT work on every possible file format
     |      supported by Bio.SeqIO (e.g. some are for multiple sequences only,
     |      and binary formats are not supported).
     |  
     |  islower(self)
     |      Return True if all ASCII characters in the record's sequence are lowercase.
     |      
     |      If there are no cased characters, the method returns False.
     |  
     |  isupper(self)
     |      Return True if all ASCII characters in the record's sequence are uppercase.
     |      
     |      If there are no cased characters, the method returns False.
     |  
     |  lower(self)
     |      Return a copy of the record with a lower case sequence.
     |      
     |      All the annotation is preserved unchanged. e.g.
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Fasta/aster.pro", "fasta")
     |      >>> print(record.format("fasta"))
     |      >gi|3298468|dbj|BAA31520.1| SAMIPF
     |      GGHVNPAVTFGAFVGGNITLLRGIVYIIAQLLGSTVACLLLKFVTNDMAVGVFSLSAGVG
     |      VTNALVFEIVMTFGLVYTVYATAIDPKKGSLGTIAPIAIGFIVGANI
     |      <BLANKLINE>
     |      >>> print(record.lower().format("fasta"))
     |      >gi|3298468|dbj|BAA31520.1| SAMIPF
     |      gghvnpavtfgafvggnitllrgivyiiaqllgstvaclllkfvtndmavgvfslsagvg
     |      vtnalvfeivmtfglvytvyataidpkkgslgtiapiaigfivgani
     |      <BLANKLINE>
     |      
     |      To take a more annotation rich example,
     |      
     |      >>> from Bio import SeqIO
     |      >>> old = SeqIO.read("EMBL/TRBG361.embl", "embl")
     |      >>> len(old.features)
     |      3
     |      >>> new = old.lower()
     |      >>> len(old.features) == len(new.features)
     |      True
     |      >>> old.annotations["organism"] == new.annotations["organism"]
     |      True
     |      >>> old.dbxrefs == new.dbxrefs
     |      True
     |  
     |  reverse_complement(self, id=False, name=False, description=False, features=True, annotations=False, letter_annotations=True, dbxrefs=False)
     |      Return new SeqRecord with reverse complement sequence.
     |      
     |      By default the new record does NOT preserve the sequence identifier,
     |      name, description, general annotation or database cross-references -
     |      these are unlikely to apply to the reversed sequence.
     |      
     |      You can specify the returned record's id, name and description as
     |      strings, or True to keep that of the parent, or False for a default.
     |      
     |      You can specify the returned record's features with a list of
     |      SeqFeature objects, or True to keep that of the parent, or False to
     |      omit them. The default is to keep the original features (with the
     |      strand and locations adjusted).
     |      
     |      You can also specify both the returned record's annotations and
     |      letter_annotations as dictionaries, True to keep that of the parent,
     |      or False to omit them. The default is to keep the original
     |      annotations (with the letter annotations reversed).
     |      
     |      To show what happens to the pre-letter annotations, consider an
     |      example Solexa variant FASTQ file with a single entry, which we'll
     |      read in as a SeqRecord:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
     |      >>> print("%s %s" % (record.id, record.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(record.letter_annotations))
     |      ['solexa_quality']
     |      >>> print(record.letter_annotations["solexa_quality"])
     |      [40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5]
     |      
     |      Now take the reverse complement, here we explicitly give a new
     |      identifier (the old identifier with a suffix):
     |      
     |      >>> rc_record = record.reverse_complement(id=record.id + "_rc")
     |      >>> print("%s %s" % (rc_record.id, rc_record.seq))
     |      slxa_0001_1_0001_01_rc NNNNNNACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
     |      
     |      Notice that the per-letter-annotations have also been reversed,
     |      although this may not be appropriate for all cases.
     |      
     |      >>> print(rc_record.letter_annotations["solexa_quality"])
     |      [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]
     |      
     |      Now for the features, we need a different example. Parsing a GenBank
     |      file is probably the easiest way to get an nice example with features
     |      in it...
     |      
     |      >>> from Bio import SeqIO
     |      >>> with open("GenBank/pBAD30.gb") as handle:
     |      ...     plasmid = SeqIO.read(handle, "gb")
     |      >>> print("%s %i" % (plasmid.id, len(plasmid)))
     |      pBAD30 4923
     |      >>> plasmid.seq
     |      Seq('GCTAGCGGAGTGTATACTGGCTTACTATGTTGGCACTGATGAGGGTGTCAGTGA...ATG')
     |      >>> len(plasmid.features)
     |      13
     |      
     |      Now, let's take the reverse complement of this whole plasmid:
     |      
     |      >>> rc_plasmid = plasmid.reverse_complement(id=plasmid.id+"_rc")
     |      >>> print("%s %i" % (rc_plasmid.id, len(rc_plasmid)))
     |      pBAD30_rc 4923
     |      >>> rc_plasmid.seq
     |      Seq('CATGGGCAAATATTATACGCAAGGCGACAAGGTGCTGATGCCGCTGGCGATTCA...AGC')
     |      >>> len(rc_plasmid.features)
     |      13
     |      
     |      Let's compare the first CDS feature - it has gone from being the
     |      second feature (index 1) to the second last feature (index -2), its
     |      strand has changed, and the location switched round.
     |      
     |      >>> print(plasmid.features[1])
     |      type: CDS
     |      location: [1081:1960](-)
     |      qualifiers:
     |          Key: label, Value: ['araC']
     |          Key: note, Value: ['araC regulator of the arabinose BAD promoter']
     |          Key: vntifkey, Value: ['4']
     |      <BLANKLINE>
     |      >>> print(rc_plasmid.features[-2])
     |      type: CDS
     |      location: [2963:3842](+)
     |      qualifiers:
     |          Key: label, Value: ['araC']
     |          Key: note, Value: ['araC regulator of the arabinose BAD promoter']
     |          Key: vntifkey, Value: ['4']
     |      <BLANKLINE>
     |      
     |      You can check this new location, based on the length of the plasmid:
     |      
     |      >>> len(plasmid) - 1081
     |      3842
     |      >>> len(plasmid) - 1960
     |      2963
     |      
     |      Note that if the SeqFeature annotation includes any strand specific
     |      information (e.g. base changes for a SNP), this information is not
     |      amended, and would need correction after the reverse complement.
     |      
     |      Note trying to reverse complement a protein SeqRecord raises an
     |      exception:
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> protein_rec = SeqRecord(Seq("MAIVMGR"), id="Test",
     |      ...                         annotations={"molecule_type": "protein"})
     |      >>> protein_rec.reverse_complement()
     |      Traceback (most recent call last):
     |         ...
     |      ValueError: Proteins do not have complements!
     |      
     |      If you have RNA without any U bases, it must be annotated as RNA
     |      otherwise it will be treated as DNA by default with A mapped to T:
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> rna1 = SeqRecord(Seq("ACG"), id="Test")
     |      >>> rna2 = SeqRecord(Seq("ACG"), id="Test", annotations={"molecule_type": "RNA"})
     |      >>> print(rna1.reverse_complement(id="RC", description="unk").format("fasta"))
     |      >RC unk
     |      CGT
     |      <BLANKLINE>
     |      >>> print(rna2.reverse_complement(id="RC", description="RNA").format("fasta"))
     |      >RC RNA
     |      CGU
     |      <BLANKLINE>
     |      
     |      Also note you can reverse complement a SeqRecord using a MutableSeq:
     |      
     |      >>> from Bio.Seq import MutableSeq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> rec = SeqRecord(MutableSeq("ACGT"), id="Test")
     |      >>> rec.seq[0] = "T"
     |      >>> print("%s %s" % (rec.id, rec.seq))
     |      Test TCGT
     |      >>> rc = rec.reverse_complement(id=True)
     |      >>> print("%s %s" % (rc.id, rc.seq))
     |      Test ACGA
     |  
     |  translate(self, table='Standard', stop_symbol='*', to_stop=False, cds=False, gap=None, id=False, name=False, description=False, features=False, annotations=False, letter_annotations=False, dbxrefs=False)
     |      Return new SeqRecord with translated sequence.
     |      
     |      This calls the record's .seq.translate() method (which describes
     |      the translation related arguments, like table for the genetic code),
     |      
     |      By default the new record does NOT preserve the sequence identifier,
     |      name, description, general annotation or database cross-references -
     |      these are unlikely to apply to the translated sequence.
     |      
     |      You can specify the returned record's id, name and description as
     |      strings, or True to keep that of the parent, or False for a default.
     |      
     |      You can specify the returned record's features with a list of
     |      SeqFeature objects, or False (default) to omit them.
     |      
     |      You can also specify both the returned record's annotations and
     |      letter_annotations as dictionaries, True to keep that of the parent
     |      (annotations only), or False (default) to omit them.
     |      
     |      e.g. Loading a FASTA gene and translating it,
     |      
     |      >>> from Bio import SeqIO
     |      >>> gene_record = SeqIO.read("Fasta/sweetpea.nu", "fasta")
     |      >>> print(gene_record.format("fasta"))
     |      >gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA) gene, partial cds
     |      CAGGCTGCGCGGTTTCTATTTATGAAGAACAAGGTCCGTATGATAGTTGATTGTCATGCA
     |      AAACATGTGAAGGTTCTTCAAGACGAAAAACTCCCATTTGATTTGACTCTGTGCGGTTCG
     |      ACCTTAAGAGCTCCACATAGTTGCCATTTGCAGTACATGGCTAACATGGATTCAATTGCT
     |      TCATTGGTTATGGCAGTGGTCGTCAATGACAGCGATGAAGATGGAGATAGCCGTGACGCA
     |      GTTCTACCACAAAAGAAAAAGAGACTTTGGGGTTTGGTAGTTTGTCATAACACTACTCCG
     |      AGGTTTGTT
     |      <BLANKLINE>
     |      
     |      And now translating the record, specifying the new ID and description:
     |      
     |      >>> protein_record = gene_record.translate(table=11,
     |      ...                                        id="phya",
     |      ...                                        description="translation")
     |      >>> print(protein_record.format("fasta"))
     |      >phya translation
     |      QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
     |      SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV
     |      <BLANKLINE>
     |  
     |  upper(self)
     |      Return a copy of the record with an upper case sequence.
     |      
     |      All the annotation is preserved unchanged. e.g.
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> record = SeqRecord(Seq("acgtACGT"), id="Test",
     |      ...                    description = "Made up for this example")
     |      >>> record.letter_annotations["phred_quality"] = [1, 2, 3, 4, 5, 6, 7, 8]
     |      >>> print(record.upper().format("fastq"))
     |      @Test Made up for this example
     |      ACGTACGT
     |      +
     |      "#$%&'()
     |      <BLANKLINE>
     |      
     |      Naturally, there is a matching lower method:
     |      
     |      >>> print(record.lower().format("fastq"))
     |      @Test Made up for this example
     |      acgtacgt
     |      +
     |      "#$%&'()
     |      <BLANKLINE>
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
     |  
     |  letter_annotations
     |      Dictionary of per-letter-annotation for the sequence.
     |      
     |      For example, this can hold quality scores used in FASTQ or QUAL files.
     |      Consider this example using Bio.SeqIO to read in an example Solexa
     |      variant FASTQ file as a SeqRecord:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
     |      >>> print("%s %s" % (record.id, record.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(record.letter_annotations))
     |      ['solexa_quality']
     |      >>> print(record.letter_annotations["solexa_quality"])
     |      [40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5]
     |      
     |      The letter_annotations get sliced automatically if you slice the
     |      parent SeqRecord, for example taking the last ten bases:
     |      
     |      >>> sub_record = record[-10:]
     |      >>> print("%s %s" % (sub_record.id, sub_record.seq))
     |      slxa_0001_1_0001_01 ACGTNNNNNN
     |      >>> print(sub_record.letter_annotations["solexa_quality"])
     |      [4, 3, 2, 1, 0, -1, -2, -3, -4, -5]
     |      
     |      Any python sequence (i.e. list, tuple or string) can be recorded in
     |      the SeqRecord's letter_annotations dictionary as long as the length
     |      matches that of the SeqRecord's sequence.  e.g.
     |      
     |      >>> len(sub_record.letter_annotations)
     |      1
     |      >>> sub_record.letter_annotations["dummy"] = "abcdefghij"
     |      >>> len(sub_record.letter_annotations)
     |      2
     |      
     |      You can delete entries from the letter_annotations dictionary as usual:
     |      
     |      >>> del sub_record.letter_annotations["solexa_quality"]
     |      >>> sub_record.letter_annotations
     |      {'dummy': 'abcdefghij'}
     |      
     |      You can completely clear the dictionary easily as follows:
     |      
     |      >>> sub_record.letter_annotations = {}
     |      >>> sub_record.letter_annotations
     |      {}
     |      
     |      Note that if replacing the record's sequence with a sequence of a
     |      different length you must first clear the letter_annotations dict.
     |  
     |  seq
     |      The sequence itself, as a Seq or MutableSeq object.
     |  
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |  
     |  __hash__ = None
    



```python
# create a Seq record 
from Bio.Seq import Seq
```


```python
simple_seq = Seq("GATC")
```


```python
simple_seq_r = SeqRecord(simple_seq)
```


```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
# we can add an ID to a seq
simple_seq_r.id = "AC1235"
```


```python
# we can add a description to a seq
simple_seq_r.description = "Made up sequence for the VDB Computtaional Biology Class"
```


```python
print(simple_seq_r.description)
```

    Made up sequence for the VDB Computtaional Biology Class



```python
simple_seq_r.seq
```




    Seq('GATC')




```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC1235', name='<unknown name>', description='Made up sequence for the VDB Computtaional Biology Class', dbxrefs=[])




```python
# how to add an annotation to a seq (called additional data)
simple_seq_r.annotations["evidence"] = "None. This is just an example"
```


```python
print(simple_seq_r.annotations["evidence"])
```

    None. This is just an example



```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC1235', name='<unknown name>', description='Made up sequence for the VDB Computtaional Biology Class', dbxrefs=[])




```python
# can do per letter annotations 
simple_seq_r.letter_annotations["phred_quality"] = [40, 40, 38, 30]
```


```python
print(simple_seq_r.letter_annotations)
```

    {'phred_quality': [40, 40, 38, 30]}



```python
# https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_00516.fna
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("NC_005816.fna.txt", "fasta")
```


```python
# displays everything that is on file for this seq 
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='gi|45478711|ref|NC_005816.1|', name='gi|45478711|ref|NC_005816.1|', description='gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
# displays record sequence 
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
# displays record ID 
record.id 
```




    'gi|45478711|ref|NC_005816.1|'




```python
# displays record description 
record.description
```




    'gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
# if it has a fasta file, this ill show that information 
record.dbxrefs
```




    []




```python
# shows this seq has no annotations
record.annotations
```




    {}




```python
record.features
```




    []




```python
# https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.gb
```


```python
record = SeqIO.read("NC_005816.gb.txt","genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
# displays the record seq ID 
record.id
```




    'NC_005816.1'




```python
record.name
```




    'NC_005816'




```python
record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
# github does not provide this info 
record.letter_annotations
```




    {}




```python
len(record.annotations)
```




    13




```python
record.annotations["source"]
```




    'Yersinia pestis biovar Microtus str. 91001'




```python
record.dbxrefs
```




    ['Project:58037']




```python
len(record.features)
```




    41




```python
from Bio import SeqFeature
```


```python
# can give a reference to where a nucleotide/ start codon could be 
start_pos = SeqFeature.AfterPosition(5)
```


```python
# can give a range for where a nucleotide/seq could be if you don't have concrete proof yet  
end_pos = SeqFeature.BetweenPosition(9, left =8, right =9)
```


```python
my_location = SeqFeature.SimpleLocation(start_pos, end_pos)
```


```python
# this shows it is somewhere before 5 to between 8 and 9
print(my_location)
```

    [>5:(8^9)]



```python
my_location.start
```




    AfterPosition(5)




```python
my_location.end
```




    BetweenPosition(9, left=8, right=9)




```python
# gives you the position in the seq where this particular section ends 
int(my_location.end)
```




    9




```python
# gives you the position in the seq where this particular section starts 
int(my_location.start)
```




    5




```python
# you cna pass the numbers into verbs 
exact_location = SeqFeature.SimpleLocation(5,9)
```


```python
print(exact_location)
```

    [5:9]



```python
exact_location.start
```




    ExactPosition(5)




```python
from Bio.SeqRecord import SeqRecord
```


```python
record = SeqRecord(Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGSETHLDSMVGQALFGD"
...         "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
...         "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
...         "SSAC"),
                    id="gi|1415038|gb|AAK54648.1|AF376133_1",
...    description="chalcone synthase [Cucumis sativus]",
... )
```


```python
# github does not do fasta files, so this is how we have to upload(what is above) and print them 
print(record.format("fasta"))
```

    >gi|1415038|gb|AAK54648.1|AF376133_1 chalcone synthase [Cucumis sativus]
    MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGSETHLDSMVGQALFGDG
    AGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISKN
    IEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNMS
    SAC
    



```python
# printing the seq without stating it as a fasta file 
print(record)
```

    ID: gi|1415038|gb|AAK54648.1|AF376133_1
    Name: <unknown name>
    Description: chalcone synthase [Cucumis sativus]
    Number of features: 0
    Seq('MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGSETHLDSMVGQ...SAC')



```python
from Bio import SeqIO
```


```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
# length displays how many base pairs 
len(record)
```




    9609




```python
len(record.features)
```




    41




```python
print(record.features[20])
```

    type: gene
    location: [4342:4780](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
print(record.features[21])
```

    type: CDS
    location: [4342:4780](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
# we can subdivide the record 
sub_record = record[4300:4800]
```


```python
len(sub_record)
```




    500




```python
len(sub_record.features)
```




    2




```python
# remember that python starts with 0 so if to loo kup the first feature this is used 
sub_record.features[0]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='gene', qualifiers=...)




```python
sub_record.features[1]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='CDS', qualifiers=...)




```python
# takes hat you had from the entire record and condenses the info down into what we 
# wanted to specifically look at in the sub_record for feature 1 
print(sub_record.features[0])
```

    type: gene
    location: [42:480](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
# for feature 2 sub_record 
print(sub_record.features[1])
```

    type: CDS
    location: [42:480](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
# can also annotate sub_records 
sub_record.annotations
```




    {'molecule_type': 'DNA'}




```python
sub_record.dbxrefs
```




    []




```python
# adds this to our annotations that we want/have for sub_record
sub_record.annotations["topology"] = "linear"
```


```python
sub_record.annotations
```




    {'molecule_type': 'DNA', 'topology': 'linear'}




```python
sub_record.id
```




    'NC_005816.1'




```python
sub_record.name
```




    'NC_005816'




```python
# shows incorrect description since for us it is not a complete seq 
sub_record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
# how to change the seq description 
sub_record.description = 'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial sequence'
```


```python
sub_record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial sequence'




```python
# we have the genbank record for this 
# displays: how long, what type, linear, when first sequenced or when first sample collected 
# it can be seen where we edited the definition to partial 
print(sub_record.format("genbank")[:200] + "...")
```

    LOCUS       NC_005816                500 bp    DNA     linear   UNK 01-JAN-1980
    DEFINITION  Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial
                sequence.
    ACCESSION   NC_00581...



```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
# this is a bacterial/ circular genome 
len(record)
```




    9609




```python
len(record.features)
```




    41




```python
record.dbxrefs
```




    ['Project:58037']




```python
record.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
# since it is a circle it's actually endless but we can shift that 
shifted = record[2000:] + record[:2000]
```


```python
shifted 
```




    SeqRecord(seq=Seq('GATACGCAGTCATATTTTTTACACAATTCTCTAATCCCGACAAGGTCGTAGGTC...GGA'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
# does not change naything within the seq, it just changes where we start/stop (creates a line instead of a circle to read)
len(shifted)
```




    9609




```python
# it is possible to lose a feature because of shifting 
len(shifted.features)
```




    40




```python
# caused because it doesn't want to carry on these annotations if they're not correct 
# it is rather safe than sorry by getting rid of them 
shifted.annotations.keys()
```




    dict_keys(['molecule_type'])




```python
shifted.dbxrefs
```




    []




```python
# if you want to keep the database cross-reference or the annotations dictionary 
# it must be done explicitly
shifted.dbxrefs = record.dbxrefs[:]
```


```python
shifted.dbxrefs
```




    ['Project:58037']




```python
shifted.annotations = record.annotations.copy()
```


```python
shifted.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
# reverse complement again 
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
# this code says to: print the first value as a string( that's what % means)
# the next  i's represent integers to print 
print("%s %i %i %i %i" % (record.id, len(record), len(record.features), len(record.dbxrefs), len(record.annotations)))
```

    NC_005816.1 9609 41 1 13



```python
rc = record.reverse_complement(id = "Testing")
```


```python
rc
```




    SeqRecord(seq=Seq('CAGGGGTCGGGGTACGCATTCCCTCATGCGTCAATATTATCTGGCATTGCGATG...ACA'), id='Testing', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
# we keep our length and features but no references or annotations
# this is done by default to make sure people don't create duplicate IDs with sequences already manipulated 
print("%s %i %i %i %i" % (rc.id, len(rc), len(rc.features), len(rc.dbxrefs), len(rc.annotations)))
```

    Testing 9609 41 0 0



```python

```


```python

```
