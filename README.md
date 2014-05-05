SeqAlign.jl
===========

Simple program to align sequences in Julia


Overview
--------

This code provides two main functions, `smithwaterman(seq1::String,
seq2::String`, and `needlemanwunsch(seq1::String, seq2::String)`
(respectively local sequence alignment, global sequence alignment).
Allowed `seq1,seq2` are just for DNA,RNA sequences (ACGT,acgt,ACGU,
acgu, or mixed capitalization). Score matrix is built in. So far the
program compute equivalent optimal solutions only for local alignments
(smithwaterman) but not yet for global alignment (needlemanwunsh).

The output is `Array{(String,String,String),1}` where the first two
strings represent the alignment, and the third a vector of symbols: 
`+ = match, o = not match, - = insert` 

Usage
-----

Just type `using SeqAlign` at julia prompt


Todos
-----

Protein sequences, custom score matrices.
 
