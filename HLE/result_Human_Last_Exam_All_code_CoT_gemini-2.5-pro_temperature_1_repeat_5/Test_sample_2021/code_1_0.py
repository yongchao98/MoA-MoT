# Define constants and initial parameters
self_pollination_rate = 0.05
cross_pollination_rate = 0.95

# The parent plant is heterozygous (Ee) for the resistance allele (e).
# Resistance is conferred by having at least one 'e' allele (Ee or ee).
# The general population is assumed to be homozygous wild-type (EE).

# 1. Calculation for self-pollination (Ee x Ee)
# Offspring: 1 EE (susceptible), 2 Ee (resistant), 1 ee (resistant)
# Probability of resistant offspring from selfing = (2 + 1) / 4 = 0.75
prob_resistant_selfing = 0.75

# 2. Calculation for cross-pollination (Ee x EE)
# Offspring: 1 EE (susceptible), 1 Ee (resistant)
# Probability of resistant offspring from crossing = 1 / 2 = 0.50
prob_resistant_crossing = 0.50

# 3. Total theoretical percentage of resistant offspring
total_resistant_fraction = (self_pollination_rate * prob_resistant_selfing) + (cross_pollination_rate * prob_resistant_crossing)
total_resistant_percentage = total_resistant_fraction * 100

# 4. Calculation for the mass increase of the E3ub protein
insertion_nucleotides = 105
nucleotides_per_codon = 3
# Average molar mass of an amino acid is ~110 g/mol, which corresponds to 110 Daltons (Da)
avg_amino_acid_mass_da = 110
mass_increase_da = (insertion_nucleotides / nucleotides_per_codon) * avg_amino_acid_mass_da
mass_increase_kda = mass_increase_da / 1000

# Print the analysis step-by-step
print("Step-by-step Analysis:")
print("="*25)

print("\n1. Calculating the percentage of drought-resistant offspring:")
print(f"The self-pollination rate is {self_pollination_rate*100}% and the cross-pollination rate is {cross_pollination_rate*100}%.")
print("For self-pollination (Ee x Ee), the probability of resistant offspring is 75%.")
print("For cross-pollination (Ee x EE), the probability of resistant offspring is 50%.")
print("\nFinal Equation for Resistance Percentage:")
print(f"Total Resistant % = (Self-pollination Rate * Resistance from Selfing) + (Cross-pollination Rate * Resistance from Crossing)")
print(f"Total Resistant % = ({self_pollination_rate} * {prob_resistant_selfing}) + ({cross_pollination_rate} * {prob_resistant_crossing})")
print(f"Total Resistant % = {self_pollination_rate * prob_resistant_selfing} + {cross_pollination_rate * prob_resistant_crossing}")
print(f"Total Resistant % = {total_resistant_fraction} or {total_resistant_percentage:.2f}%")

print("\n2. Analyzing protein activity:")
print("Co-expression of Par22 with E3ub-wt decreased Par22 levels, indicating E3ub-wt is an active E3 ubiquitin ligase.")
print("Co-expression with E3ub-insert105 increased Par22 levels, indicating the insertion inactivates the protein's ligase function.")
print("Conclusion: Only E3ub-wt is an active ubiquitin ligase.")

print("\n3. Analyzing protein-protein interaction:")
print("Native mass spectrometry of Par22 (50kDa) + E3ub-wt (60kDa) showed a complex at 110kDa, confirming interaction.")
print("The mixture of Par22 + E3ub-insert105 showed separate peaks for each protein, indicating no interaction.")
print("Conclusion: Par22 cannot interact with E3ub-insert105.")


print("\n4. Calculating the protein mass increase from the insertion:")
print(f"The insertion has {insertion_nucleotides} nucleotides.")
print("\nFinal Equation for Mass Increase:")
print(f"Mass Increase (kDa) = (Insertion Nucleotides / Nucleotides per Codon) * Avg. Amino Acid Mass (Da) / 1000")
print(f"Mass Increase (kDa) = ({insertion_nucleotides} / {nucleotides_per_codon}) * {avg_amino_acid_mass_da} / 1000")
print(f"Mass Increase (kDa) = {insertion_nucleotides / nucleotides_per_codon:.1f} * {avg_amino_acid_mass_da} / 1000")
print(f"Mass Increase (kDa) = {mass_increase_da} / 1000 = {mass_increase_kda:.2f} kDa")
print("This is approximately 4.0 kDa.")


print("\n" + "="*25)
print("Final Summary:")
print(f"- Theoretically, {total_resistant_percentage:.2f}% of the offspring should be drought-resistant.")
print("- Only E3ub-wt is an active ubiquitin ligase.")
print("- Par22 cannot interact with E3ub-insert105.")
print(f"- The insertion increases the mass of E3ub by approximately 4.0 kDa.")
print("\nThis combination of findings corresponds to option J.")