import math

# Step 1: Calculate the mass increase of the E3ub protein.
nucleotide_insertion_length = 105
nucleotides_per_codon = 3
avg_mass_amino_acid_da = 110  # Daltons
da_per_kda = 1000

# The number of amino acids added by the insertion.
added_amino_acids = nucleotide_insertion_length / nucleotides_per_codon

# The mass increase in kDa.
mass_increase_da = added_amino_acids * avg_mass_amino_acid_da
mass_increase_kda = mass_increase_da / da_per_kda

# Step 2 & 3: Analyze experimental data (this is a reasoning step).
# - Activity: The wild-type E3ub-wt reduces Par22 levels (from 700 to 200), so it is an active ubiquitin ligase.
#   The E3ub-insert105 variant does not (levels increase to 3000), so it is inactive.
# - Interaction: E3ub-wt (60kDa) and Par22 (50kDa) form a 110kDa complex, so they interact.
#   E3ub-insert105 and Par22 show separate peaks (50kDa for Par22), meaning they do not interact.

# Step 4: Calculate the theoretical percentage of resistant offspring.
self_pollination_rate = 0.05
cross_pollination_rate = 0.95

# For self-pollination (Wt/Ins x Wt/Ins), offspring are 1/4 Wt/Wt, 2/4 Wt/Ins, 1/4 Ins/Ins.
# Assuming the 'Ins' allele is dominant, Wt/Ins and Ins/Ins are resistant.
# Probability of resistant offspring from selfing is 3/4.
prob_resistant_selfing = 3 / 4

# For cross-pollination with the general susceptible population (Wt/Ins x Wt/Wt),
# offspring are 1/2 Wt/Wt, 1/2 Wt/Ins.
# Probability of resistant offspring from crossing is 1/2.
prob_resistant_crossing = 1 / 2

# Total probability of resistant offspring is the weighted average.
total_resistant_prob = (self_pollination_rate * prob_resistant_selfing) + (cross_pollination_rate * prob_resistant_crossing)
total_resistant_percent = total_resistant_prob * 100

# Step 5: Print the final conclusion, showing the numbers used in the equation.
print("Summary of Findings:")
print(f"1. The 105 nucleotide insertion increases the protein mass by approximately {mass_increase_kda:.1f} kDa.")
print("2. Experimental data shows that only E3ub-wt is an active ubiquitin ligase.")
print("3. Experimental data also shows that Par22 cannot interact with E3ub-insert105.")
print("\nCalculation of Resistant Offspring Percentage:")
print(f"The total percentage is calculated as: ({self_pollination_rate} * {prob_resistant_selfing}) + ({cross_pollination_rate} * {prob_resistant_crossing}) = {total_resistant_prob}")
print(f"This gives a total of {total_resistant_percent:.2f}% resistant offspring.")

print("\nFinal Conclusion:")
print(f"Theoretically, {total_resistant_percent:.2f}% of the offspring should be drought-resistant. "
      "Only E3ub-wt is an active ubiquitin ligase. "
      "Par22 cannot interact with E3ub-insert105. "
      f"The insertion increases the mass of E3ub by approximately {mass_increase_kda:.1f} kDA.")