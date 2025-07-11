# Step 1: Calculate the theoretical offspring resistance rate.
self_pollination_rate = 0.05
cross_pollination_rate = 0.95

# For self-pollination (Ee x Ee), resistant offspring (EE, Ee) are 3/4.
resistance_from_selfing = 0.75
# For cross-pollination (Ee x ee), resistant offspring (Ee) are 1/2.
resistance_from_crossing = 0.50

# Total theoretical resistance rate
total_resistance = (self_pollination_rate * resistance_from_selfing) + (cross_pollination_rate * resistance_from_crossing)
total_resistance_percent = total_resistance * 100

print("--- Offspring Resistance Calculation ---")
print(f"Contribution from self-pollination (5% of cases): {self_pollination_rate} * {resistance_from_selfing} = {self_pollination_rate * resistance_from_selfing}")
print(f"Contribution from cross-pollination (95% of cases): {cross_pollination_rate} * {resistance_from_crossing} = {cross_pollination_rate * resistance_from_crossing}")
print(f"Total theoretical resistance rate: {total_resistance} or {total_resistance_percent:.2f}%\n")


# Step 2 & 3: Analyze protein activity and interaction.
print("--- Protein Function Analysis ---")
print("Co-expression results show that E3ub-wt reduces Par22 levels (700 -> 200 units), so it is an active E3 ligase.")
print("Co-expression with E3ub-insert105 does not reduce Par22 levels (700 -> 3000 units), so it is NOT an active E3 ligase.")
print("Native MS shows Par22 (50kDa) and E3ub-wt (60kDa) form a 110kDa complex, so they interact.")
print("Native MS with E3ub-insert105 shows a peak for Par22 alone (50kDa), so they do NOT interact.\n")


# Step 4: Calculate the mass increase from the insertion.
insertion_length_nt = 105
codons_per_amino_acid = 3
avg_mass_amino_acid_Da = 110

# Calculate the number of added amino acids
added_amino_acids = insertion_length_nt / codons_per_amino_acid

# Calculate the mass increase in Daltons and kiloDaltons
mass_increase_Da = added_amino_acids * avg_mass_amino_acid_Da
mass_increase_kDa = mass_increase_Da / 1000

print("--- Mass Increase Calculation ---")
print(f"The insertion of {insertion_length_nt} nucleotides adds {int(added_amino_acids)} amino acids.")
print(f"The estimated mass increase is {int(added_amino_acids)} * {avg_mass_amino_acid_Da} Da = {mass_increase_Da:.0f} Da.")
print(f"This is approximately {mass_increase_kDa:.1f} kDa.\n")

# Step 5: Synthesize the final conclusion matching the correct option.
print("--- Final Conclusion ---")
print(f"Theoretically, {total_resistance_percent:.2f}% of the offspring should be drought-resistant.")
print("Only E3ub-wt is an active ubiquitin ligase.")
print("Par22 cannot interact with E3ub-insert105.")
print(f"The insertion increases the mass of E3ub by approximately {mass_increase_kDa:.1f}kDA.")
print("This matches option J.")

<<<J>>>