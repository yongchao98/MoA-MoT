import math

# Step 1: Calculate the theoretical percentage of resistant offspring.
# Define the given probabilities and ratios.
self_pollination_rate = 0.05
cross_pollination_rate = 0.95

# In self-pollination (wt/ins x wt/ins), resistant genotypes are wt/ins and ins/ins.
# This corresponds to 3/4 of the offspring.
resistant_ratio_selfing = 0.75

# In cross-pollination (wt/ins x wt/wt), the only resistant genotype is wt/ins.
# This corresponds to 1/2 of the offspring.
resistant_ratio_crossing = 0.50

# Calculate the total proportion of resistant offspring.
total_resistant_proportion = (self_pollination_rate * resistant_ratio_selfing) + (cross_pollination_rate * resistant_ratio_crossing)
total_resistant_percentage = total_resistant_proportion * 100

print("Step 1: Calculate the percentage of resistant offspring.")
print(f"The total percentage of resistant offspring is calculated as: ({self_pollination_rate} * {resistant_ratio_selfing}) + ({cross_pollination_rate} * {resistant_ratio_crossing}) = {total_resistant_proportion:.4f}")
print(f"This equals {total_resistant_percentage:.2f}%.\n")

# Step 2: Analyze the E3ub ligase activity from the experimental data.
print("Step 2: Analyze E3ub ligase activity.")
print("The amount of Par22 drops from 700 to 200 units when co-expressed with E3ub-wt, indicating Par22 degradation. Therefore, E3ub-wt is an active ubiquitin ligase.")
print("The amount of Par22 increases from 700 to 3000 units when co-expressed with E3ub-insert105, indicating no degradation. Therefore, E3ub-insert105 is inactive.\n")

# Step 3: Analyze the protein interaction from mass spectrometry data.
print("Step 3: Analyze protein interaction.")
print("Mass spectrometry of Par22 (50kDa) and E3ub-wt (60kDa) shows a peak at 110kDa, which is the sum of their masses. This confirms they interact.")
print("Mass spectrometry of Par22 and E3ub-insert105 shows a peak at 50kDa, corresponding to free Par22. This confirms they do not interact.\n")

# Step 4: Calculate the mass increase from the nucleotide insertion.
nucleotide_insertion_length = 105
codons_per_amino_acid = 3
avg_amino_acid_mass_da = 110 # in Daltons

# Calculate the number of amino acids added.
added_amino_acids = nucleotide_insertion_length / codons_per_amino_acid

# Calculate the total mass increase in Daltons.
mass_increase_da = added_amino_acids * avg_amino_acid_mass_da

# Convert the mass increase to kiloDaltons (kDa).
mass_increase_kda = mass_increase_da / 1000

print("Step 4: Calculate the mass increase of the E3ub protein.")
print(f"The insertion has {nucleotide_insertion_length} nucleotides, which code for {int(added_amino_acids)} amino acids ({nucleotide_insertion_length} / {codons_per_amino_acid}).")
print(f"The mass increase is approximately {int(added_amino_acids)} amino acids * {avg_amino_acid_mass_da} Da/amino acid = {int(mass_increase_da)} Da.")
print(f"This is equal to {mass_increase_kda:.2f} kDa, which is approximately 4.0 kDa.\n")

# Step 5: Final Conclusion
print("Conclusion:")
print("Based on the analysis, the correct statements are:")
print(f"- Theoretically, {total_resistant_percentage:.2f}% of the offspring should be drought-resistant.")
print("- Only E3ub-wt is an active ubiquitin ligase.")
print("- Par22 cannot interact with E3ub-insert105.")
print(f"- The insertion increases the mass of E3ub by approximately {round(mass_increase_kda, 1)} kDA.")
print("This matches all claims in option J.")
