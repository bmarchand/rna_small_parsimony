# rna_small_parsimony
Some code to experiment around the different notions of distances to solve the small parsimony problem on RNA structures.

# Steps to reproduce results

```
pip install .
```

To install the package rnadist, which includes useful routines.

Then, with the Rfam.seed file placed in the experiments/rf_small_parsimony/resources folder, execute in that same folder 

```
python3 extract_seed_alignments.py
```

Also, retrieve from RFAM the ``tree files'' zip and extract it into the tree_files subfolder.

Finally, in the experiments/rf_small_parsimony folder, execute

```
snakemake -c4
```

To produce all the data, and

```
python3 average_num_bps.py
```

To produce the figure.
