# NT-AcPredictor

NT-AcPredictor: predicting N-terminal acetylated sequences by decision tree

## Abstract

This program predicts whether N-terminal acetylation occurs for proteins whose initiator methionine residues are removed by methionine amino peptidase.

The program was produced and would be maintained by Kazunori D Yamada (kyamada@ecei.tohoku.ac.jp) and Masaru Miyagi (mxm356@case.edu).


## Usage

To run the program, simply type the following command. The format of input file must be FASTA format. The program accepts multiple sequence file.

$ nT-AcPredictor.py -i INPUT -f (truncated|exist)


## Reference

Yamada KD, Omori S, Nish H, Miyagi M, Identification of the sequence determinants of protein N-terminal acetylation through a decision tree approach, BMC Bioinformatics, 18(1):289, 2017
