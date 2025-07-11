# 1. Define variables from the problem description
self_pollination_rate = 0.05
cross_pollination_rate = 0.95 # This is 1.0 - self_pollination_rate

# Resistance probability for each pollination type
# Self-pollination (Wt/Ins x Wt/Ins) -> 1/4 Wt/Wt, 2/4 Wt/Ins, 1/4 Ins/Ins
# Resistant offspring are Wt/Ins and Ins/Ins, so 3/4 are resistant.
resistance_from_selfing = 0.75

# Cross-pollination (Wt/Ins x Wt/Wt) -> 1/2 Wt/Wt, 1/2 Wt/Ins
# Resistant offspring are Wt/Ins, so 1/2 are resistant.
resistance_from_crossing = 0.50

# 2. Calculate the total percentage of resistant offspring
prob_resistant = (self_pollination_rate * resistance_from_selfing) + (cross_pollination_rate * resistance_from_crossing)
percent_resistant = prob_resistant * 100

print("Calculation for resistant offspring percentage:")
print(f"({self_pollination_rate} * {resistance_from_selfing}) + ({cross_pollination_rate} * {resistance_from_crossing}) = {prob_resistant}")
print(f"Total percentage of resistant offspring: {percent_resistant:.2f}%\n")


# 3. Calculate the protein mass increase from the nucleotide insertion
insertion_length_nt = 105
nucleotides_per_amino_acid = 3
avg_mass_per_amino_acid_da = 110 # in Daltons

# Calculate number of added amino acids
num_added_amino_acids = insertion_length_nt / nucleotides_per_amino_acid

# Calculate mass increase in Daltons and kiloDaltons
mass_increase_da = num_added_amino_acids * avg_mass_per_amino_acid_da
mass_increase_kda = mass_increase_da / 1000

print("Calculation for protein mass increase:")
print(f"{insertion_length_nt} nucleotides / {nucleotides_per_amino_acid} nt_per_aa = {int(num_added_amino_acids)} amino acids")
print(f"{int(num_added_amino_acids)} amino acids * {avg_mass_per_amino_acid_da} Da/aa = {int(mass_increase_da)} Da")
print(f"Mass increase: {mass_increase_kda:.2f} kDa, which is approximately 4.0 kDa.\n")

# 4. Synthesize the final conclusion based on calculations and experimental data
print("Final Conclusion:")
print(f"Theoretically, {percent_resistant:.2f}% of the offspring should be drought-resistant. "
      "Only E3ub-wt is an active ubiquitin ligase. "
      "Par22 cannot interact with E3ub-insert105. "
      f"The insertion increases the mass of E3ub by approximately {mass_increase_kda:.1f}kDA.")