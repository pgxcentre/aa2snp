# Introduction

Queries the Ensembl REST API to convert an amino acid change to a SNP. I did
not do a lot of checking (_e.g._ the script could return synonymous variants if
they overlap with the protein position).
**Use at your own risk**

To use the script, clone the repository and run it using Python 3. I only use 
modules from the standard library so this should not be a problem.

Running the script with ``--help`` should provide enough information to use it.
Here is an example:

```console
$ python3 aa2snp.py -g ENSP00000369497 -v F101L
id	chromosome	position	strand	alleles	amino_changes
rs397507301	13	32893449	1	C/A	F/L
```

