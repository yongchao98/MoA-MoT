# Part 1: Calculate the percentage of drought-resistant offspring
self_pollination_rate = 0.05
cross_pollination_rate = 0.95
# In a self-cross (Ee x Ee), resistant offspring (Ee, ee) are 3/4 of the total.
resistant_proportion_selfing = 0.75
# In a cross with wild-type (Ee x EE), resistant offspring (Ee) are 1/2 of the total.
resistant_proportion_crossing = 0.50

# Part 2: Calculate protein mass increase from the insertion
insertion_nucleotide_length = 105
# 3 nucleotides code for 1 amino acid
added_amino_acids = insertion_nucleotide_length / 3
# Average mass of an amino acid is ~110 Daltons
avg_amino_acid_mass_da = 110
mass_increase_da = added_amino_acids * avg_amino_acid_mass_da
# Convert Daltons to kiloDaltons
mass_increase_kda = mass_increase_da / 1000

print("--- Analysis Calculations ---")

# Print Offspring Percentage Calculation
print("\n1. Calculation of Drought-Resistant Offspring Percentage:")
total_resistant_fraction = (self_pollination_rate * resistant_proportion_selfing) + (cross_pollination_rate * resistant_proportion_crossing)
print(f"The equation for the total fraction of resistant offspring is: ({self_pollination_rate} * {resistant_proportion_selfing}) + ({cross_pollination_rate} * {resistant_proportion_crossing})")
print(f"This evaluates to: {total_resistant_fraction:.4f}, or {total_resistant_fraction * 100:.2f}%")

# Print Mass Increase Calculation
print("\n2. Calculation of Protein Mass Increase:")
print(f"Number of added amino acids = {insertion_nucleotide_length} nucleotides / 3 = {int(added_amino_acids)} amino acids")
print(f"Mass increase in Daltons = {int(added_amino_acids)} amino acids * {avg_amino_acid_mass_da} Da/amino acid = {int(mass_increase_da)} Da")
print(f"Mass increase in kiloDaltons = {int(mass_increase_da)} Da / 1000 = {mass_increase_kda:.2f} kDa")

print("\n--- Final Conclusion based on Analysis ---")
print("The calculated results and interpretation of the experimental data align perfectly with answer choice J.")