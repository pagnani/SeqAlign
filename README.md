SeqAlign.jl
===========

Simple program to align sequences in Julia


Overview
--------

This code provides two main functions, `smithwaterman(seq1::String,
seq2::String`, and `needlemanwunsch(seq1::String, seq2::String)`
(respectively local sequence alignment, global sequence alignment).
Allowed `seq1,seq2` are just for DNA,RNA sequences (ACGT,acgt,ACUT,
acut). Score matrix is built in.


Usage
-----

Just type `using SeqAlign` at julia prompt


Todos
-----

Protein sequences, custom score matrices.
 
