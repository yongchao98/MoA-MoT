# Step 1: Calculate the protein mass increase
nucleotide_insertion_length = 105
codons_per_amino_acid = 3
avg_amino_acid_mass_da = 110

# Calculate the number of amino acids added
added_amino_acids = nucleotide_insertion_length / codons_per_amino_acid

# Calculate the total mass increase in Daltons and convert to kiloDaltons (kDa)
mass_increase_da = added_amino_acids * avg_amino_acid_mass_da
mass_increase_kda = mass_increase_da / 1000

print(f"Analysis of Insertion Mass:")
print(f"The insertion of {nucleotide_insertion_length} nucleotides adds {int(added_amino_acids)} amino acids.")
print(f"The estimated protein mass increase is {mass_increase_da} Da, or approximately {mass_increase_kda:.1f} kDa.")
print("-" * 20)

# Step 2 & 3 are based on interpretation of the experimental data provided in the text.
# Interpretation of Activity (Co-expression experiment):
# - Par22 alone: 700 units (baseline)
# - Par22 + E3ub-wt: 200 units (level decreased -> E3ub-wt is active)
# - Par22 + E3ub-insert105: 3000 units (level increased -> E3ub-insert105 is inactive)
# Conclusion: Only E3ub-wt is an active ubiquitin ligase.

# Interpretation of Interaction (Mass Spectrometry):
# - Par22 (50kDa) + E3ub-wt (60kDa) -> 110kDa peak (interaction occurs)
# - Par22 (50kDa) + E3ub-insert105 -> 50kDa peak (free Par22) and another peak (no complex peak) -> No interaction occurs.
# Conclusion: Par22 cannot interact with E3ub-insert105.

# Step 4: Calculate the theoretical percentage of resistant offspring
self_pollination_rate = 0.05
cross_pollination_rate = 0.95

# In self-pollination (W/I x W/I), offspring genotypes are 1/4 W/W, 2/4 W/I, 1/4 I/I.
# Resistant offspring have at least one I allele (W/I or I/I).
resistance_rate_selfing = 2/4 + 1/4

# In cross-pollination (W/I x W/W), offspring are 1/2 W/W, 1/2 W/I.
# Resistant offspring are W/I.
resistance_rate_crossing = 1/2

# Calculate the total resistance frequency
total_resistance = (resistance_rate_selfing * self_pollination_rate) + (resistance_rate_crossing * cross_pollination_rate)

# Convert to percentage
total_resistance_percent = total_resistance * 100

print("Calculation of Resistant Offspring Percentage:")
print(f"The self-pollination rate is {self_pollination_rate * 100}%. The probability of resistant offspring from selfing is {resistance_rate_selfing}.")
print(f"The cross-pollination rate is {cross_pollination_rate * 100}%. The probability of resistant offspring from crossing is {resistance_rate_crossing}.")
print(f"Final Equation: ({resistance_rate_selfing} * {self_pollination_rate}) + ({resistance_rate_crossing} * {cross_pollination_rate}) = {total_resistance}")
print(f"Theoretically, {total_resistance_percent:.2f}% of the offspring should be drought-resistant.")
print("-" * 20)

print("Summary of conclusions:")
print("1. Mass Increase: ~4.0 kDa")
print("2. Activity: Only E3ub-wt is an active ubiquitin ligase.")
print("3. Interaction: Par22 cannot interact with E3ub-insert105.")
print("4. Offspring Resistance: 51.25%")
print("\nMatching these conclusions with the options points to option J.")
