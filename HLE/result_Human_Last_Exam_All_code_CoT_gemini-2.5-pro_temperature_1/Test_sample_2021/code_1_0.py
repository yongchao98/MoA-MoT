import math

# --- Introduction ---
print("This script analyzes the provided biological data to find the correct conclusion.")
print("It performs two main calculations: the protein mass increase and the theoretical percentage of resistant offspring.")
print("-" * 40)

# --- Part 1: Calculate the Mass Increase ---
print("Step 1: Calculating the mass increase from the 105 nucleotide insertion.")

inserted_nucleotides = 105
nucleotides_per_codon = 3
# The average molecular weight of an amino acid residue in a protein is ~110 Daltons.
avg_mass_amino_acid_Da = 110

# Calculate the number of added amino acids
added_amino_acids = inserted_nucleotides / nucleotides_per_codon

# Calculate the total mass increase in Daltons and kiloDaltons
mass_increase_Da = added_amino_acids * avg_mass_amino_acid_Da
mass_increase_kDa = mass_increase_Da / 1000

print(f"The insertion of {inserted_nucleotides} nucleotides adds {int(added_amino_acids)} amino acids.")
print(f"The calculated mass increase is {mass_increase_Da:.0f} Da, which is approximately {round(mass_increase_kDa, 1)} kDa (~4.0 kDa).")
print("-" * 40)

# --- Part 2: Calculate Theoretical Offspring Resistance ---
print("Step 2: Calculating the theoretical percentage of drought-resistant offspring.")

prob_self_pollination = 0.05
prob_cross_pollination = 0.95

# Resistance fraction from self-pollination of a heterozygote (e.g., Rr x Rr)
# Offspring are 1/4 RR, 2/4 Rr, 1/4 rr. Resistant genotypes are RR and Rr.
resistance_from_selfing = 3 / 4

# Resistance fraction from cross-pollination with a homozygous recessive (e.g., Rr x rr)
# We assume the general population is wild-type (non-resistant).
# Offspring are 1/2 Rr, 1/2 rr. Resistant genotype is Rr.
resistance_from_crossing = 1 / 2

# Calculate the final percentage using the law of total probability
total_resistant_fraction = (prob_self_pollination * resistance_from_selfing) + (prob_cross_pollination * resistance_from_crossing)
total_resistant_percentage = total_resistant_fraction * 100

print("The calculation for the total percentage of resistant offspring is a weighted average based on pollination rates:")
print(f"Final Equation: ({prob_self_pollination} * {resistance_from_selfing}) + ({prob_cross_pollination} * {resistance_from_crossing}) = {total_resistant_fraction}")
print(f"This results in a theoretical resistant population of {total_resistant_percentage:.2f}%.")
print("-" * 40)

# --- Part 3: Conclusion ---
print("Summary of Findings:")
print(f"1. Offspring Resistance: {total_resistant_percentage:.2f}%")
print("2. Ligase Activity: Only E3ub-wt is an active ubiquitin ligase (based on Par22 degradation).")
print("3. Protein Interaction: Only E3ub-wt interacts with Par22 (based on mass spectrometry complex formation).")
print(f"4. Mass Increase: The insertion adds approximately {round(mass_increase_kDa, 1)} kDa.")
print("\nThese findings align with option J.")
