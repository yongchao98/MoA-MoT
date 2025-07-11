import sys

# Step 1: Calculate the theoretical percentage of drought-resistant offspring.

# Define pollination rates
self_pollination_rate = 0.05
cross_pollination_rate = 0.95

# Define resistance rates for each type of pollination
# For self-pollination (WI x WI -> 1/4 WW, 2/4 WI, 1/4 II), 3/4 or 75% are resistant (WI, II)
resistance_from_self = 0.75
# For cross-pollination (WI x WW -> 1/2 WW, 1/2 WI), 1/2 or 50% are resistant (WI)
resistance_from_cross = 0.50

# Calculate the total proportion of resistant offspring
total_resistance = (self_pollination_rate * resistance_from_self) + (cross_pollination_rate * resistance_from_cross)
resistance_percentage = total_resistance * 100

print("--- Analysis of Offspring Resistance ---")
print(f"The parent plant is heterozygous (WI). Resistance is conferred by the 'I' allele.")
print(f"Contribution from self-pollination ({self_pollination_rate*100}%): {self_pollination_rate} * {resistance_from_self} = {self_pollination_rate * resistance_from_self}")
print(f"Contribution from cross-pollination ({cross_pollination_rate*100}%): {cross_pollination_rate} * {resistance_from_cross} = {cross_pollination_rate * resistance_from_cross}")
print(f"Total theoretical resistance = ({self_pollination_rate} * {resistance_from_self}) + ({cross_pollination_rate} * {resistance_from_cross}) = {total_resistance}")
print(f"Therefore, theoretically {resistance_percentage:.2f}% of the offspring should be drought-resistant.\n")

# Step 2: Analyze protein function (E3 ligase activity).
print("--- Analysis of Protein Function ---")
print("In the co-expression experiment, Par22 levels dropped from 700 to 200 units with E3ub-wt, indicating degradation.")
print("This means E3ub-wt is an active ubiquitin ligase.")
print("With E3ub-insert105, Par22 levels increased to 3000 units, indicating no degradation.")
print("This means E3ub-insert105 is NOT an active ubiquitin ligase.\n")

# Step 3: Analyze protein-protein interaction.
print("--- Analysis of Protein Interaction ---")
print("In native mass spectrometry, Par22 (50kDa) and E3ub-wt (60kDa) showed a single peak at 110kDa.")
print("Since 50 + 60 = 110, this confirms a direct interaction.")
print("Par22 (50kDa) and E3ub-insert105 showed separate peaks, not a combined one.")
print("This confirms that Par22 cannot interact with E3ub-insert105.\n")


# Step 4: Calculate the mass increase from the insertion.
insertion_length_nt = 105
codons_per_aa = 3
avg_mass_aa_da = 115 # Average mass of an amino acid in Daltons

# Calculate the number of amino acids added
num_added_aa = insertion_length_nt / codons_per_aa

# Calculate the mass increase
mass_increase_da = num_added_aa * avg_mass_aa_da
mass_increase_kda = mass_increase_da / 1000

print("--- Analysis of Mass Increase ---")
print(f"The insertion is {insertion_length_nt} nucleotides long.")
print(f"Number of amino acids added = {insertion_length_nt} / {codons_per_aa} = {int(num_added_aa)}")
print(f"Estimated mass increase = {int(num_added_aa)} amino acids * {avg_mass_aa_da} Da/aa = {mass_increase_da:.0f} Da")
print(f"This is approximately {mass_increase_kda:.1f} kDa.\n")

# Step 5: Final Conclusion
print("--- Final Conclusion ---")
print(f"Based on the analysis:")
print(f"- Theoretically, {resistance_percentage:.2f}% of the offspring should be drought-resistant.")
print(f"- Only E3ub-wt is an active ubiquitin ligase.")
print(f"- Par22 cannot interact with E3ub-insert105.")
print(f"- The insertion increases the mass of E3ub by approximately {mass_increase_kda:.1f}kDA.")
print("\nThis matches option J.")
sys.stdout = open('/dev/null', 'w') # Suppress further output
#<<<J>>>