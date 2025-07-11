# The sequence for the intracellular linker domain between TM3 and TM4
# of the human GABAArho1 receptor (UniProt P24046, residues 325-427).
tm3_tm4_linker_sequence = "LFDRLVNCVFENPFKNDLMYWVIETAFSLPTAGFSSRLALHNKQKARALEGLRASPALQLWVKKAIPDLVFMFISILCVYLAFTYHVHVSERKPERAAKARI"

# Count the number of Cysteine ('C') residues in the sequence.
cysteine_count = tm3_tm4_linker_sequence.count('C')

# Find the positions (1-based index) of each Cysteine residue in the full protein.
# The linker starts at position 325.
start_position = 325
cysteine_positions = []
for i, amino_acid in enumerate(tm3_tm4_linker_sequence):
    if amino_acid == 'C':
        # Add the position relative to the full protein sequence.
        cysteine_positions.append(start_position + i)

# Print the results
print(f"The TM3-TM4 linker sequence is: {tm3_tm4_linker_sequence}")
print(f"The Cysteine residues are found at the following positions in the full protein: {cysteine_positions}")
print(f"The total number of Cysteine residues in the TM3-TM4 linker is: {cysteine_count}")