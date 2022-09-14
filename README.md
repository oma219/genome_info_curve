# Investigating New Genomic Information in Each New Genome

There are various measures that are used to estimate the amount of distinct information in a text. In particular, the measures that we investigate are the r (number
of runs in Burrows-Wheeler Transform of text) and delta (maximum ratio of number of distinct kmers over k).

Some things `generate_curve.py` does for us.

* Sweeping k's for kmc
* Intermediate file wrangling
* Adding reverse complement

We're choosing to leave plasmids in for now.

### Prep

```
git clone https://github.com/oma219/genome_info_curve
git clone https://github.com/alshai/pfbwt-f
cd pfbwt-f
git submodule init
git submodule update
mamba install seqtk
```

Check that you can run:

```
/scratch16/blangme2/langmead/PFP_LZ77/build/_deps/bigbwt-build/newscanNT.x
/scratch16/blangme2/langmead/PFP_LZ77/build/src/lz_77_test
```

### To compute revcomps

```
seqtk seq -r <FASTA_IN> | sed 's/^>/>revcomp_/' > <FASTA_OUT_TMP>
cat <FASTA_IN> <FASTA_OUT_TMP> > <FASTA_OUT>
```

### To get r

```
./pfbwt-f/build/pfbwt-f64 -o <OUT_PREFIX> <FASTA>
```

r is printed in the output

### To get delta

Per k:

```
kmc -ci1 -cs1 -k30 -b -t12 -fm <FASTA> <OUT_PREFIX> <WORKING_DIR>
```

number of unique k-mers is printed in the output

### To get z

Just use default PFP parameters

```
~/scr16_blangme2/langmead/PFP_LZ77/build/_deps/bigbwt-build/newscanNT.x <FASTA> -w 10 -p 100 -f
~/scr16_blangme2/langmead/PFP_LZ77/build/src/lz_77_test <FASTA>
```

z is printed in the output under "phrase number:  XYZ"

### Input format

Perhaps need to preprocess all the files to both exclude the header and add the reverse complement.

Actually, we think they all understand FASTA headers because they use kseq.  So maybe all we need to do is add a record with the reverse complement before indexing.

https://github.com/lh3/seqtk/blob/master/kseq.h
