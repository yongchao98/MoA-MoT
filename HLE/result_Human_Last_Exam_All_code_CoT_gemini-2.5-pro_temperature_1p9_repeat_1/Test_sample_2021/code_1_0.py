# Step 1: Analyze function (logic based on problem description)
# - Co-expression: Par22 level drops with E3ub-wt (degradation) but not with E3ub-insert105.
#   Conclusion: Only E3ub-wt is an active ubiquitin ligase.
# - Mass Spec: Par22 (50kDa) + E3ub-wt (60kDa) -> 110kDa complex.
#   Conclusion: E3ub-wt interacts with Par22.
# - Mass Spec: Par22 (50kDa) + E3ub-insert105 -> 50kDa peak (free Par22).
#   Conclusion: E3ub-insert105 does not interact with Par22.
# Overall Conclusion 1: Only E3ub-wt is active, and only E3ub-wt interacts with Par22.

# Step 2: Calculate the mass increase from the nucleotide insertion
insertion_length_bp = 105
codons_per_amino_acid = 3
avg_amino_acid_mass_da = 110  # Daltons

# Calculate the number of inserted amino acids
num_inserted_amino_acids = insertion_length_bp / codons_per_amino_acid

# Calculate the total mass increase in Daltons
mass_increase_da = num_inserted_amino_acids * avg_amino_acid_mass_da

# Convert mass to kiloDaltons (kDa)
mass_increase_kda = mass_increase_da / 1000

print("--- Mass Increase Calculation ---")
print(f"The insertion of {insertion_length_bp} bp corresponds to {int(num_inserted_amino_acids)} amino acids.")
print(f"The estimated mass increase is {int(num_inserted_amino_acids)} aa * {avg_amino_acid_mass_da} Da/aa = {int(mass_increase_da)} Da.")
print(f"This is approximately {mass_increase_kda:.2f} kDa, which rounds to 4.0 kDa.\n")

# Step 3: Calculate the percentage of resistant offspring
self_pollination_rate = 0.05
cross_pollination_rate = 0.95

# For self-pollination (Rw x Rw), offspring are 1/4 RR, 1/2 Rw, 1/4 ww.
# Resistant offspring (RR + Rw) are 3/4 or 75%.
prob_resistant_selfing = 0.75

# For cross-pollination (Rw x ww), offspring are 1/2 Rw, 1/2 ww.
# Resistant offspring (Rw) are 1/2 or 50%.
prob_resistant_crossing = 0.50

# Calculate the contribution from each pollination type
contribution_selfing = self_pollination_rate * prob_resistant_selfing
contribution_crossing = cross_pollination_rate * prob_resistant_crossing

# Calculate the total proportion of resistant offspring
total_resistant_proportion = contribution_selfing + contribution_crossing
total_resistant_percentage = total_resistant_proportion * 100

print("--- Offspring Resistance Calculation ---")
print("Equation for total resistant offspring proportion:")
print(f"({self_pollination_rate} * {prob_resistant_selfing}) + ({cross_pollination_rate} * {prob_resistant_crossing}) = {total_resistant_proportion}")
print(f"The theoretical percentage of drought-resistant offspring is {total_resistant_percentage:.2f}%\n")

# Step 4: Synthesize conclusions
print("--- Final Conclusion ---")
print(f"1. Offspring Resistance: {total_resistant_percentage:.2f}%")
print("2. Enzyme Activity: Only E3ub-wt is an active ubiquitin ligase.")
print("3. Protein Interaction: Only E3ub-wt can interact with Par22.")
print(f"4. Mass Increase: The insertion adds approximately {round(mass_increase_kda,1)} kDa to the protein mass.")
print("\nThese results align with option A.")