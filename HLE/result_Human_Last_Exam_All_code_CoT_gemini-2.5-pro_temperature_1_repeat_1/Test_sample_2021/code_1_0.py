import math

# Part 1: Calculate the theoretical percentage of resistant offspring.
self_pollination_rate = 0.05
cross_pollination_rate = 0.95

# For self-pollination (Ee x Ee), resistant offspring (Ee, ee) are 3/4 of the progeny.
resistance_from_selfing = 3/4

# For cross-pollination (Ee x EE), resistant offspring (Ee) are 1/2 of the progeny.
resistance_from_crossing = 1/2

# The final equation for the total resistance frequency:
# (self_pollination_rate * resistance_from_selfing) + (cross_pollination_rate * resistance_from_crossing)
total_resistance_fraction = (self_pollination_rate * resistance_from_selfing) + (cross_pollination_rate * resistance_from_crossing)
total_resistance_percentage = total_resistance_fraction * 100

# Part 2: Analyze protein function and interaction from the text.
# E3ub-wt decreases Par22 levels, so it is an active ligase.
# E3ub-insert105 increases Par22 levels, so it is inactive.
conclusion_activity = "Only E3ub-wt is an active ubiquitin ligase"

# Native MS shows a complex for E3ub-wt + Par22.
# Native MS shows free Par22 for E3ub-insert105 + Par22, meaning no stable interaction.
conclusion_interaction = "Par22 cannot interact with E3ub-insert105"

# Part 3: Calculate the mass increase of the protein.
insertion_length_bp = 105
codons_per_aa = 3
avg_aa_mass_da = 110 # Daltons

# The final equation for the mass increase:
# (insertion_length_bp / codons_per_aa) * avg_aa_mass_da
num_amino_acids = insertion_length_bp / codons_per_aa
mass_increase_da = num_amino_acids * avg_aa_mass_da
mass_increase_kda = mass_increase_da / 1000
# The calculated mass is 3.85 kDa, which is approximately 4.0 kDa as stated in the answer choice.
mass_increase_kda_approx = 4.0

# Part 4: Assemble and print the final conclusion, which corresponds to option J.
print("Based on the data analysis:")
# Print each number in the final equations as requested
print(f"The calculation for offspring resistance is: ({self_pollination_rate} * {resistance_from_selfing}) + ({cross_pollination_rate} * {resistance_from_crossing}) = {total_resistance_fraction:.4f} or {total_resistance_percentage:.2f}%")
print(f"The calculation for protein mass increase is: ({insertion_length_bp} / {codons_per_aa}) * {avg_aa_mass_da} Da = {mass_increase_kda:.2f} kDa, which is approximately {mass_increase_kda_approx:.1f} kDa.")
print("\nFinal Conclusion:")
print(f"Theoretically, {total_resistance_percentage:.2f}% of the offspring should be drought-resistant. {conclusion_activity}. {conclusion_interaction}. The insertion increases the mass of E3ub by approximately {mass_increase_kda_approx:.1f}kDA.")
<<<J>>>