# --- Define constants based on the problem description ---
self_pollination_rate = 0.05
cross_pollination_rate = 0.95
insertion_sequence = "gatgaaaaagtgtttaccaaagaactggatcagaccattgaacagctgaacgaatgcggccagctgagcgaaagccaggtgaaaagcctgatgtgcgaagcgacc"
avg_amino_acid_mass_da = 110

# --- 1. Calculate the percentage of resistant offspring ---
print("Step 1: Calculating the theoretical percentage of resistant offspring.")
# For self-pollination (WI x WI), resistant genotypes are WI and II.
# The ratio is 1 WW : 2 WI : 1 II. So, 3 out of 4 (75%) are resistant.
resistant_fraction_selfing = 0.75
# For cross-pollination (WI x WW), resistant genotype is WI.
# The ratio is 1 WW : 1 WI. So, 1 out of 2 (50%) are resistant.
resistant_fraction_crossing = 0.50

# Calculate the weighted average
total_resistant_percentage = (resistant_fraction_selfing * self_pollination_rate + resistant_fraction_crossing * cross_pollination_rate) * 100

print("The calculation for the percentage of resistant offspring is:")
print(f"({resistant_fraction_selfing} * {self_pollination_rate}) + ({resistant_fraction_crossing} * {cross_pollination_rate}) = {total_resistant_percentage / 100}")
print(f"Therefore, theoretically {total_resistant_percentage:.2f}% of the offspring should be drought-resistant.\n")

# --- 2. Analyze the protein data and calculate the mass increase ---
print("Step 2: Analyzing protein function and calculating the mass increase from the insertion.")
# Logic based on experimental data
print("From the co-expression data:")
print("- E3ub-wt reduces Par22 levels (700 -> 200), indicating it is an active E3 ubiquitin ligase that targets Par22 for degradation.")
print("- E3ub-insert105 increases Par22 levels (700 -> 3000), indicating it is inactive.\n")

print("From the native mass spectrometry data:")
print("- E3ub-wt (60kDa) and Par22 (50kDa) form a 110kDa complex, showing they interact.")
print("- With E3ub-insert105, a peak for free Par22 (50kDa) is observed, showing Par22 cannot interact with E3ub-insert105.\n")

# Calculate mass increase
insertion_length_bp = len(insertion_sequence)
added_amino_acids = insertion_length_bp / 3
mass_increase_da = added_amino_acids * avg_amino_acid_mass_da
mass_increase_kda = mass_increase_da / 1000

print("The calculation for the mass increase of the E3ub protein is:")
print(f"Insertion length = {insertion_length_bp} nucleotides")
print(f"Number of added amino acids = {insertion_length_bp} / 3 = {added_amino_acids:.0f}")
print(f"Mass increase (Da) = {added_amino_acids:.0f} amino acids * {avg_amino_acid_mass_da} Da/amino acid = {mass_increase_da:.0f} Da")
print(f"Mass increase (kDa) = {mass_increase_da:.0f} Da / 1000 = {mass_increase_kda:.2f} kDa")
print(f"The insertion increases the mass of E3ub by approximately {mass_increase_kda:.1f} kDa.\n")

print("Summary of conclusions:")
print(f"- Offspring resistance: {total_resistant_percentage:.2f}%")
print("- E3ub activity: Only E3ub-wt is active.")
print("- Interaction: Par22 cannot interact with E3ub-insert105.")
print(f"- Mass increase: ~{mass_increase_kda:.1f} kDa.")
print("\nThese conclusions match answer choice J.")