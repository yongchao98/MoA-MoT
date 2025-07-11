import math

# Step 1: Define variables from the problem statement
self_pollination_rate = 0.05
cross_pollination_rate = 0.95
insertion_sequence = "gatgaaaaagtgtttaccaaagaactggatcagaccattgaacagctgaacgaatgcggccagctgagcgaaagccaggtgaaaagcctgatgtgcgaagcgacc"
avg_aa_mass_da = 110 # Average mass of an amino acid in Daltons

# Step 2: Analyze Protein Function and Interaction based on experimental data
# This part is based on interpretation of the text, not calculation.
# - E3ub-wt activity: Degrades Par22 (200 < 700), so it is an active ligase.
# - E3ub-insert105 activity: Stabilizes Par22 (3000 > 700), so it is inactive.
# - E3ub-wt interaction: Forms a 110kDa complex (50kDa + 60kDa), so it interacts with Par22.
# - E3ub-insert105 interaction: Results in separate protein peaks, so it does NOT interact with Par22.

# Step 3: Calculate the mass increase of the E3ub protein
insertion_length_nt = len(insertion_sequence)
codons_added = insertion_length_nt / 3
mass_increase_da = codons_added * avg_aa_mass_da
mass_increase_kda = mass_increase_da / 1000

# Step 4: Calculate the theoretical percentage of resistant offspring
# The parent plant is heterozygous (Genotype: Rr, where 'R' is the resistance allele).
# Resistance is conferred by the R allele, so Rr and RR genotypes are resistant.

# For self-pollination (Rr x Rr), offspring are 1 RR : 2 Rr : 1 rr.
# Resistant fraction = (RR + Rr) / Total = 3/4
resistance_from_selfing = 3/4

# For cross-pollination with the general population (assumed wild-type, rr),
# the cross is Rr x rr. Offspring are 1 Rr : 1 rr.
# Resistant fraction = Rr / Total = 1/2
resistance_from_crossing = 1/2

# Total resistant offspring proportion is a weighted average
total_resistant_proportion = (self_pollination_rate * resistance_from_selfing) + (cross_pollination_rate * resistance_from_crossing)
total_resistant_percentage = total_resistant_proportion * 100

# Step 5: Print the full analysis and conclusion
print("Analysis Summary:")
print("1. Protein Activity: Based on the co-expression data, only E3ub-wt is an active ubiquitin ligase.")
print("2. Protein Interaction: Based on mass spectrometry, only E3ub-wt interacts with Par22. Par22 cannot interact with E3ub-insert105.")
print("-" * 40)

print("Calculation of Protein Mass Increase:")
print(f"The insertion of {insertion_length_nt} nucleotides adds {int(codons_added)} amino acids.")
print(f"Calculated mass increase: {int(codons_added)} aa * {avg_aa_mass_da} Da/aa = {mass_increase_da} Da, or {mass_increase_kda:.1f} kDa (approx. 4.0 kDa).")
print("-" * 40)

print("Calculation of Resistant Offspring Percentage:")
print("The final equation for the proportion of resistant offspring is:")
print(f"({self_pollination_rate} * {resistance_from_selfing}) + ({cross_pollination_rate} * {resistance_from_crossing}) = {total_resistant_proportion}")
print(f"The theoretical percentage of drought-resistant offspring is {total_resistant_percentage:.2f}%.")
print("-" * 40)

print("\nFinal Conclusion:")
print("The results show that 51.25% of the offspring should be drought-resistant, only E3ub-wt is an active ligase, Par22 cannot interact with E3ub-insert105, and the insertion increases E3ub mass by approximately 4.0 kDa. This aligns perfectly with choice J.")

<<<J>>>