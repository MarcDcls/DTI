Author : Marc Duclusaud \
Contact : marc.duclusaud@u-bordeaux.fr 

# Deaminase Target Identifier (DTI)

## Description

This package is a tool to identify possible targets of a deaminase in a set of sequences of nucleotides. The tool identify target codons in proximity with PAMs (Protospacer Adjacent Motifs) in the sequences. The tool also checks the presence of off-targets in a genome file.

## Installation

Installing the module tqdm allows to vizualise the execution in the terminal, while the module multiprocessing allows to speed up the identification process by using multiple CPU cores. To install these modules, run the following command:

```bash
pip install tqdm multiprocessing
```


## Usage

Create a folder in the data repository and place the fasta and genome files in it. The arborescence should look like this:

```
data/
  example/
    example.fasta  # fasta file containing the sequences of nucleotides to be identified
    example.fa     # genome file to check the presence of off-targets
```

Then update the identification parameters in the file `identifier.py`. The parameters are:
- **targets** : target codons for the deaminase with the index of the target nucleotide in the codon
- **pams** : PAMs (Protospacer Adjacent Motifs)
- **deaminase_window** : window for the nucleotide targetted by the deaminase before the PAM
- **concordance_threshold** : threshold of concordance for a sequence to be considered as an off-target
- **allowed_mismatches** : number of mismatches allowed in an off-target
- **allowed_gaps** : number of gaps allowed in an off-target
- **nb_processes** : number of processes to use for the identification (default is the number of CPU cores available)

Finally run the following command:

```bash
python identifier.py
```

## Targets

Targets are defined as a list containing the following elements:
- The target codon (3 nucleotides) to be identified
- The index of the target nucleotide in the codon (0, 1 or 2)
- The strand on which the deaminase will be used (coding or complementary)

**Example 1:** ["CAG", 0, "coding"] means that the target is the 1st nucleotide (C) of the codon CAG. To target this nucleotide, the cytosine deaminase will be used on the coding strand.

**Example 2:** ["TGA", 1, "complementary"] means that the target is the 2nd nucleotide (G) of the codon TGA. To target this nucleotide,  the cytosine deaminase need to target a C nucleotide, and thus will be used on the complementary strand.

## Output

The output is a text result file created in the folder containing the fasta and genome files. It contains for each entry of the fasta file the sequences of nucleotides identified as possible targets, as well as the potential off-targets found in the genome file.

**Note:** In an off-target, the mismatches are indicated by "*" and the gaps by "-".

## Example

The repository `data/CDS_mbov_PG45_example` contains a fasta file `CDS_mbov_PG45_example.fasta` and a genome file `Mbov_PG45.fa`. 

To perform the identification, simply run the command: 

```bash
python identifier.py
```