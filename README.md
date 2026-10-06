# PSYLAB_ANOPLURA

## Contents & Description

### 1. Assembled Genomes

All genomes assembled in this study are archived in EMBL Flat File format (`*.embl`).  
GFF3 format annotation file (`*.gff`) and FASTA format genomic DNA (`*.fasta`) or Protein sequences (`*.pep`) also provided.  

1. S1844378 (*Pediculus schaeffi*)
2. S1844379 (*Pthirus gorillae*)
3. S1844380 (*Pthirus pubis*)

### 2. Script
1. AAAPSA (All-Against-All Pairwise Sequence Alignment)
   Multiprocessing integrated All-against-all Pairwise Alignment tool.  
   It requires Python and its module `psa`.
2. `assembly_pipeline.sh`
   Automated script to pre-processing of raw data and assembly with `SPAdes` (via `Shovil`).  
   It requires BBTools and Lighter, Shovil.
3. `repeat_masking_cmd.sh`
   For repeat masking. This script is collection of command lines orignally authored by Daren Card.
4. `fetch_busco-odb10_description.sh`
   For collecting BUSCO catalog data. This script is collection of command lines orignally authored by Daren Card.
5. `collect_complete-nonredundant_busco.sh`
   This script is collection of command lines orignally authored by Daren Card.
6. `train_augustus_by_busco.sh`
   For running AUGUSTUS `eTraining.pl`. This script is collection of command lines orignally authored by Daren Card.
