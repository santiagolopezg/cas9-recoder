#########################

File needed:

- sequence.fasta (nucleotide sequence to be processed, in fasta format)

#########################

Basic how-to:

- Clone or fork the `cas9-recoder` repo;
- Open a terminal window / python shell;
- Run `python mk_up_dst_region.py`:
	This command will generate a sequence.txt file. Containing the original fasta sequence, padded with A's so as to put the start codon at position 100.
- Run `python aminoacid_tiling.py` or `python codon_change.py`, depending on whether the desired outcome is Alanine tiling or systematic recoding of each codon.
- Run `python filter_for_pos.py`:
	This will filter all generated guide+donor pairs, and keep only the ones which contain recodings within the desired range.
- run `python sanity_check.py`:
	This checks all files, filtering the generated sequences against guides in donor sequences, or presence of PAM sequences in donors.
- run `python validation.py`:
	This performs a final twofold validation of the sequences by:
	1. aligning the donor with the reference sequence;
	2. translating it and align it with reference (also translated).
