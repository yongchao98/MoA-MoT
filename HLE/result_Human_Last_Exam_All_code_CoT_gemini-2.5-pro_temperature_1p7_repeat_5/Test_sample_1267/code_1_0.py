import random

# --- Analysis based on the problem ---
# Original anticodon 5'-UAA-3' reads mRNA codon 5'-UUA-3' (Leucine)
# Mutated anticodon 5'-UUG-3' reads mRNA codon 5'-CAA-3' (Glutamine)
# The mutated tRNA is still charged with Leucine but now reads a Glutamine codon.
# This results in competition at the CAA codon between the correct tRNA-Gln and the incorrect tRNA-Leu.

# --- Parameters from the problem statement ---
codon_in_question = "CAA"
correct_amino_acid = "Glutamine"
misincorporated_amino_acid = "Leucine"
total_instances = 1000
misincorporated_instances = 1

# --- Calculation ---
# This calculates the number of correct events based on the total and misincorporated events.
correct_instances = total_instances - misincorporated_instances
misincorporation_frequency = misincorporated_instances / total_instances

# --- Output ---
print(f"Analysis of the tRNA Mutation:")
print(f"The mutated tRNA causes a misreading of the '{codon_in_question}' mRNA codon.")
print(f"Normally, this codon specifies the amino acid '{correct_amino_acid}'.")
print(f"However, the mutated tRNA causes '{misincorporated_amino_acid}' to be inserted instead.")
print("\nThis substitution happens at a low frequency due to competition with the correct tRNA.")

print(f"\nBased on the provided rate of 1 in {total_instances}:")
# The problem asks to output each number in the final equation.
print(f"Misincorporation Frequency = {misincorporated_instances} / {total_instances} = {misincorporation_frequency}")
print(f"This means for every {total_instances} times the '{codon_in_question}' codon is translated:")
print(f"- Approximately {correct_instances} times, the correct amino acid ({correct_amino_acid}) is inserted.")
print(f"- Approximately {misincorporated_instances} time, the incorrect amino acid ({misincorporated_amino_acid}) is inserted.")
