# -Trio-Angiosperm-Sequence-Analysis-Toolkit
This script automates angiosperm trio analysis. It aligns sequences (MAFFT), calculates Tajima's D &amp; M1/M2 values, and converts to PHYLIP. Supports batch &amp; single file processing.
Angiosperm_Calculations.py
This Python script facilitates statistical analyses on aligned FASTA file datasets containing trios
of annual, perennial, and outgroup plant species. It offers functionalities for sequence alignment,
Tajima's D test execution, and PHYLIP format conversion.
Functionalities:
1. Alignment: Aligns FASTA files in batch mode (entire directory) or compiled mode (single
file).
2. Tajima's D test: Performs Tajima's D test on FASTA files in batch mode or compiled
mode. Outputs results (filename, M1, M2, D1, D2, chi-squared) to a CSV file.
3. PHYLIP conversion: Converts FASTA files to PHYLIP format (batch mode or compiled
mode) for compatibility with PAML servers.
a. Note: This functionality currently supports FASTA files where the entire sequence
for each species is on the same line directly after the '>' tag.

Command-line Arguments:
● -a, --annual: Fasta tag specifying the annual species.
● -p, --perennial: Fasta tag specifying the perennial species.
● -r, --outgroup: Fasta tag specifying the outgroup species.
● -t, --type: Type of sequences being analyzed (compiled or batch).
● -o, --output: Output directory or file (depends on the function).
● -i, --input: Input directory or file (depends on the function).
● -f, --function: Function to perform (align, tajima, phylip).
Command-line Examples:
1. Batch Alignment:
Bash
python3 Angiosperm_Calculations.py -f align -i input_dir -o output_dir -t batch

This command aligns all FASTA files within the input_dir directory and saves the aligned
outputs in the output_dir directory (batch mode).
2. Compiled Alignment:
Bash
python3 Angiosperm_Calculations.py -f align -i input_file.fasta -o path/to/output/filename(no
.phy) -t compiled

This command aligns a single FASTA file named input_file.fasta and saves the aligned
output in the output_dir directory (compiled mode).
3. Batch Tajima's D test:
Bash
python3 Angiosperm_Calculations.py -f tajima -i input_dir -o output.csv -t batch -a annual_tag -p
perennial_tag -r outgroup_tag

This command performs Tajima's D test on all FASTA files within the input_dir directory. It
outputs the results (M1, M2, D1, D2, chi-squared values) to a CSV file named output.csv.
Replace annual_tag, perennial_tag, and outgroup_tag with the corresponding Fasta
tags for your species.
4. Compiled Tajima's D test:
Bash
python3 Angiosperm_Calculations.py -f tajima -i input_file.fasta -o output.csv -t compiled -a
annual_tag -p perennial_tag -r outgroup_tag

This command performs Tajima's D test on a single FASTA file named input_file.fasta. It
outputs the results (M1, M2, D1, D2, chi-squared values) to a CSV file named output.csv.
Replace annual_tag, perennial_tag, and outgroup_tag with the corresponding Fasta
tags for your species.
5. Batch PHYLIP conversion:
Bash
python3 Angiosperm_Calculations.py -f phylip -i input_dir -o output_dir -t batch

This command converts all FASTA files within the input_dir directory to PHYLIP format and
saves them in the output_dir directory (batch mode).

6. Compiled PHYLIP conversion:
Bash
python3 Angiosperm_Calculations.py -f phylip -i input_file.fasta -o output_dir -t compiled

This command converts a single FASTA file named input_file.fasta to PHYLIP format and
saves the output file in the output_dir directory (compiled mode).
