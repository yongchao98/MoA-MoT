# Step 1: Calculate the theoretical percentage of resistant offspring.
self_pollination_rate = 0.05
cross_pollination_rate = 0.95

# The resistance allele is dominant, so for self-pollination (Wt/Ins x Wt/Ins),
# the resistant genotypes (Wt/Ins and Ins/Ins) make up 75% of the offspring.
resistant_from_selfing_fraction = 0.75

# For cross-pollination with a wild-type plant (Wt/Ins x Wt/Wt),
# the resistant genotype (Wt/Ins) makes up 50% of the offspring.
resistant_from_crossing_fraction = 0.50

# The equation for the total percentage of resistant offspring is:
# (self-pollination rate * fraction resistant from selfing) + (cross-pollination rate * fraction resistant from crossing)
total_resistant_offspring_fraction = (self_pollination_rate * resistant_from_selfing_fraction) + (cross_pollination_rate * resistant_from_crossing_fraction)
total_resistant_offspring_percentage = total_resistant_offspring_fraction * 100

print("--- Offspring Resistance Calculation ---")
print(f"Contribution from self-pollination: {self_pollination_rate} * {resistant_from_selfing_fraction} = {self_pollination_rate * resistant_from_selfing_fraction}")
print(f"Contribution from cross-pollination: {cross_pollination_rate} * {resistant_from_crossing_fraction} = {cross_pollination_rate * resistant_from_crossing_fraction}")
print(f"Total theoretical resistant offspring: ({self_pollination_rate} * {resistant_from_selfing_fraction}) + ({cross_pollination_rate} * {resistant_from_crossing_fraction}) = {total_resistant_offspring_fraction:.4f}")
print(f"This corresponds to {total_resistant_offspring_percentage:.2f}% of the offspring.\n")


# Step 2: Calculate the approximate mass increase from the nucleotide insertion.
nucleotide_insertion_length = 105
codons_per_amino_acid = 3
avg_mass_per_amino_acid_Da = 110 # Approximate average mass in Daltons

# Equation for number of amino acids added:
amino_acids_added = nucleotide_insertion_length / codons_per_amino_acid

# Equation for mass increase in Daltons:
mass_increase_Da = amino_acids_added * avg_mass_per_amino_acid_Da
mass_increase_kDa = mass_increase_Da / 1000

print("--- Protein Mass Increase Calculation ---")
print(f"Number of amino acids added: {nucleotide_insertion_length} nucleotides / {codons_per_amino_acid} = {int(amino_acids_added)} amino acids")
print(f"Approximate mass increase: {int(amino_acids_added)} amino acids * {avg_mass_per_amino_acid_Da} Da/amino acid = {int(mass_increase_Da)} Da")
print(f"This is approximately {mass_increase_kDa:.1f} kDa.\n")

# Step 3: Summarize experimental findings and select the final answer.
print("--- Summary of Findings ---")
print("1. Offspring Resistance: Calculated to be 51.25%.")
print("2. E3 Ubiquitin Ligase Activity: Only E3ub-wt reduces Par22 levels, so only E3ub-wt is active.")
print("3. Protein Interaction: Mass spectrometry shows a complex for E3ub-wt and unbound Par22 for E3ub-insert105, meaning Par22 cannot interact with E3ub-insert105.")
print("4. Mass Increase: Calculated to be approximately 4.0 kDa.")

print("\n--- Final Answer ---")
print("The analysis corresponds to option J:")
print("Theoretically, 51.25% of the offspring should be drought-resistant. Only E3ub-wt is an active ubiquitin ligase. Par22 cannot interact with E3ub-insert105. The insertion increases the mass of E3ub by approximately 4.0kDA.")
