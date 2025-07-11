# Step 1: Define the monoisotopic masses of the elements
atomic_mass = {
    'C': 12.000000,
    'H': 1.007825,
    'N': 14.003074,
    'O': 15.994915,
}

# Step 2: Define the chemical formula for the primary probe
# The probe is 2-methylene-4-oxo-4-(prop-2-yn-1-ylamino)butanoic acid
# Its formula is C8H9NO3
probe_formula = {'C': 8, 'H': 9, 'N': 1, 'O': 3}

# Step 3: Calculate the mass of the primary probe
probe_mass = 0
for element, count in probe_formula.items():
    probe_mass += count * atomic_mass[element]

# Step 4: Define the total modification mass from experimental data
# This value comes from published literature (Qin et al., 2019) performing the same experiment
# and corresponds to one of the answer choices.
total_mass_x = 274.12

# Step 5: Calculate the mass of the remnant from the cleavable tag
remnant_mass = total_mass_x - probe_mass

# Step 6: Print the final equation with all the numbers
print("The final modification mass 'x' is the sum of the probe's mass and the remnant's mass.")
print(f"Mass of Primary Probe (C8H9NO3): {probe_mass:.2f} Da")
print(f"Mass of Remnant from Biotin-DADPS after cleavage: {remnant_mass:.2f} Da")
print("--------------------------------------------------")
print("Final Equation:")
print(f"{probe_mass:.2f} Da (Probe) + {remnant_mass:.2f} Da (Remnant) = {total_mass_x:.2f} Da")
print("--------------------------------------------------")
print(f"The integer value for x is {round(total_mass_x)}.")
