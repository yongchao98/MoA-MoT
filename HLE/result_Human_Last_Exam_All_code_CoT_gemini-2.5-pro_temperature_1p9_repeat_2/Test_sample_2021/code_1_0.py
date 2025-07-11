import math

# Step 1: Calculate the percentage of drought-resistant offspring

# Given probabilities
self_pollination_rate = 0.05
cross_pollination_rate = 0.95

# The parent plant genotype is heterozygous (wt/insert).
# Assumption: The insertion allele ('insert') is dominant, so both (wt/insert) and (insert/insert) genotypes are resistant.

# Case 1: Self-pollination (wt/insert x wt/insert)
# Offspring probabilities: 0.25 wt/wt, 0.50 wt/insert, 0.25 insert/insert
# Proportion of resistant offspring in this case:
p_resistant_self = 0.50 + 0.25  # (wt/insert) + (insert/insert)

# Case 2: Cross-pollination with wild-type population (wt/insert x wt/wt)
# Offspring probabilities: 0.50 wt/wt, 0.50 wt/insert
# Proportion of resistant offspring in this case:
p_resistant_cross = 0.50 # (wt/insert)

# Total percentage of resistant offspring
total_resistant_percent = (self_pollination_rate * p_resistant_self + cross_pollination_rate * p_resistant_cross) * 100

# Step 2: Analyze protein activity and interaction from the text
# Activity Analysis:
# - Par22 + E3ub-wt: Par22 level drops from 700 to 200. This indicates degradation.
#   Conclusion: E3ub-wt is an active ubiquitin ligase.
# - Par22 + E3ub-insert105: Par22 level rises from 700 to 3000. No degradation occurs.
#   Conclusion: E3ub-insert105 is NOT an active ubiquitin ligase.
activity_conclusion = "Only E3ub-wt is an active ubiquitin ligase."

# Interaction Analysis (Mass Spectrometry):
# - Par22 (50kDa) + E3ub-wt (60kDa) -> peak at 110kDa.
#   Conclusion: The proteins form a 1:1 complex, so they interact.
# - Par22 (50kDa) + E3ub-insert105 -> peaks at 50kDa and 690kDa.
#   Conclusion: The presence of a 50kDa peak means free Par22 is present, so it does NOT interact with E3ub-insert105.
interaction_conclusion = "Par22 cannot interact with E3ub-insert105."


# Step 3: Calculate the mass increase of the protein E3ub-insert105

# Given values
insertion_length_bp = 105
bp_per_amino_acid = 3
# Average molar mass of an amino acid residue is ~110 Daltons (g/mol)
avg_amino_acid_mass_da = 110

# Calculation
inserted_amino_acids = insertion_length_bp / bp_per_amino_acid
mass_increase_da = inserted_amino_acids * avg_amino_acid_mass_da
mass_increase_kda = mass_increase_da / 1000

# Step 4: Synthesize the final answer statement

print("Based on the analysis:")
print(f"Theoretically, {total_resistant_percent:.2f}% of the offspring should be drought-resistant.")
print(f"Based on the co-expression experiment, {activity_conclusion}")
print(f"Based on the mass spectrometry data, {interaction_conclusion}")
print(f"The {insertion_length_bp} nucleotide insertion increases the mass of E3ub by approximately {mass_increase_kda:.1f} kDa.")

print("\nFinal conclusion sentence:")
print(f"Theoretically, {total_resistant_percent:.2f}% of the offspring should be drought-resistant. {activity_conclusion}. {interaction_conclusion}. The insertion increases the mass of E3ub by approximately {mass_increase_kda:.1f} kDA.")