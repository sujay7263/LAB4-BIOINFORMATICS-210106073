# LAB4-BIOINFORMATICS-210106073


**Progressive Multiple Sequence Alignment**

This Python software uses a fundamental dynamic programming technique to conduct progressive multiple-sequence alignment. It uses protein sequences supplied by the user as input and outputs the aligned sequences.
For this code to run we need to install biopython package, we can use the command "pip install biopython".

Using dynamic programming, the script implements a condensed version of the progressive multiple-sequence alignment process. There are two main steps in it:

Pair-wise Alignment: The script uses the Needleman-Wunsch method to pair-wise align sequences while taking match, mismatch, and gap penalties into account.

Profile Alignment: A profile alignment is constructed using the previously aligned sequences. The final aligned sequences are created by aligning this profile with each succeeding sequence.

