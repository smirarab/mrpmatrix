mrpmatrix
=========

This java program builds a mrpmatrix on a huge input files quickly. 

It uses some tricks to use very little time and memory.


```
Usage: <treesfile> <output> <ouputformat> [-dna] [-randomize seed]
                <treesfile>: A file containing newick trees, one tree per line
                <Output>: The name of the output MRP Matrix file
                <outformat>: use NEXUS for nexus, PHYLIP for phylip, or FASTA for fasta fromatted otuput
                -dna: output As and Ts instead of 0 and 1
                -randomize: randomize 0-1 codings. Seed number is optional.
```

The randomize option is useful if you are doing an MRL analyses, but is not needed for MRP. 
