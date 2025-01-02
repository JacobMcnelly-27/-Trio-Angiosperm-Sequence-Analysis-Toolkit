# -Trio-Angiosperm-Sequence-Analysis-Toolkit
Functionalities:

Sequence Alignment: Aligns FASTA files in batch or compiled mode (single file) using MAFFT (separate installation required).
Tajima's D Test: Performs Tajima's D test on FASTA files (batch/compiled mode) and outputs results (M1, M2, D1, D2, chi-squared) to a CSV file.
PHYLIP Conversion: Converts FASTA files to PHYLIP format (batch/compiled mode) for compatibility with PAML servers. (Note: Currently supports FASTA files with sequences on a single line after the '>' tag.)
Command-Line Arguments:

-a, --annual: Fasta tag for the annual species.
-p, --perennial: Fasta tag for the perennial species.
-r, --outgroup: Fasta tag for the outgroup species.
-t, --type: Type of sequences (compiled or batch).
-o, --output: Output directory or file (depends on function).
-i, --input: Input directory or file (depends on function).
-f, --function: Function to perform (align, tajima, phylip).
Command-Line Examples:

Batch Alignment:
Bash

python3 Angiosperm_Calculations.py -f align -i input_dir -o output_dir -t batch
Compiled Alignment:
Bash

python3 Angiosperm_Calculations.py -f align -i input_file.fasta -o path/to/output/filename(no .phy) -t compiled
Batch Tajima's D test:
Bash

python3 Angiosperm_Calculations.py -f tajima -i input_dir -o output.csv -t batch -a annual_tag -p perennial_tag -r outgroup_tag
Compiled Tajima's D test:
Bash

python3 Angiosperm_Calculations.py -f tajima -i input_file.fasta -o output.csv -t compiled -a annual_tag -p perennial_tag -r outgroup_tag
Batch PHYLIP conversion:
Bash

python3 Angiosperm_Calculations.py -f phylip -i input_dir -o output_dir -t batch
Compiled PHYLIP conversion:
Bash

python3 Angiosperm_Calculations.py -f phylip -i input_file.fasta -o output_dir -t compiled
Note: Replace annual_tag, perennial_tag, and outgroup_tag with the corresponding Fasta tags for your species.
