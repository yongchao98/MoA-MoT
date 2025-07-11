# Part 1: Calculate the mass increase from the nucleotide insertion.
nucleotide_insertion_length = 105
# Every 3 nucleotides code for 1 amino acid.
amino_acids_added = nucleotide_insertion_length / 3
# The average molecular weight of an amino acid is ~110 Daltons (Da).
average_aa_mass_da = 110
mass_increase_da = amino_acids_added * average_aa_mass_da
# Convert Daltons to kiloDaltons (kDa).
mass_increase_kda = mass_increase_da / 1000

print("Step 1: Analysis of Mass Increase")
print(f"The insertion adds {int(amino_acids_added)} amino acids.")
print(f"Calculated mass increase: {mass_increase_kda:.1f} kDa.\n")

print("Step 2: Analysis of Protein Function and Interaction")
print("E3 Ligase Activity: The sharp decrease in Par22 with E3ub-wt indicates it's an active E3 ligase. The increase with E3ub-insert105 shows it's inactive.")
print("Protein Interaction: The 110kDa complex from E3ub-wt + Par22 confirms interaction. The presence of free 50kDa Par22 with E3ub-insert105 confirms no interaction.\n")

# Part 2: Calculate the theoretical percentage of resistant offspring.
self_pollination_rate = 0.05
cross_pollination_rate = 0.95

# For self-pollination (WI x WI), resistant offspring (WI, II) are 3 out of 4.
resistant_proportion_selfing = 0.75
# For cross-pollination (WI x WW), resistant offspring (WI) are 2 out of 4.
resistant_proportion_crossing = 0.50

# Calculate the total theoretical percentage.
total_resistant_proportion = (self_pollination_rate * resistant_proportion_selfing) + (cross_pollination_rate * resistant_proportion_crossing)
total_resistant_percentage = total_resistant_proportion * 100

print("Step 3: Calculation of Resistant Offspring Percentage")
print("The final calculation for the proportion of resistant offspring is:")
# The prompt requires printing each number in the final equation.
print(f"({self_pollination_rate} * {resistant_proportion_selfing}) + ({cross_pollination_rate} * {resistant_proportion_crossing}) = {total_resistant_proportion}")
print(f"\nThis corresponds to a theoretical resistance rate of {total_resistant_percentage:.2f}%.\n")

print("Conclusion:")
print("Combining all points:")
print("- Offspring resistance is 51.25%.")
print("- Only E3ub-wt is an active ubiquitin ligase.")
print("- Par22 cannot interact with E3ub-insert105.")
print("- The insertion increases mass by approximately 4.0 kDa.")
print("This matches option J.")
