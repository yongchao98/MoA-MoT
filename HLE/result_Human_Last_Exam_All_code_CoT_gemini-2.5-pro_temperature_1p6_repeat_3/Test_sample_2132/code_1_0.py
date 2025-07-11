# Step 1: Identify the product and its properties.
# The Favorskii rearrangement of 2-chlorocyclohexanone results in a ring contraction,
# yielding the product cyclopentanecarboxylic acid.
# The chemical formula is C6H10O2.

# Step 2: Define a function to calculate Böttcher Molecular Complexity (BMC).
def calculate_bmc(name, num_atoms, num_bonds, num_rings, num_rotatable_bonds):
    """Calculates and prints the Böttcher Molecular Complexity."""
    print(f"Calculating Böttcher Molecular Complexity for: {name}")
    print(f"Molecular formula: C6H10O2\n")

    print(f"1. Number of atoms (N_atoms): {num_atoms}")
    print(f"2. Number of bonds (N_bonds): {num_bonds}")
    print(f"3. Number of rings (N_rings): {num_rings}")
    print(f"4. Number of rotatable bonds (N_rotatable_bonds): {num_rotatable_bonds}")
    
    # Check for division by zero
    if num_rotatable_bonds == 0:
        print("\nCannot calculate BMC as the number of rotatable bonds is zero.")
        return

    # Step 3: Calculate the BMC using the formula.
    # BMC = (N_atoms * N_bonds * N_rings) / N_rotatable_bonds
    bmc_value = (num_atoms * num_bonds * num_rings) / num_rotatable_bonds

    # Step 4: Print the final equation and the result.
    print("\nBMC Formula: (N_atoms * N_bonds * N_rings) / N_rotatable_bonds")
    print("\nFinal Calculation:")
    print(f"{num_atoms} * {num_bonds} * {num_rings} / {num_rotatable_bonds} = {bmc_value}")

# Parameters for cyclopentanecarboxylic acid
product_name = "Cyclopentanecarboxylic acid"

# N_atoms = 6 (Carbon) + 10 (Hydrogen) + 2 (Oxygen)
n_atoms = 18

# N_bonds = For a molecule with one ring, N_bonds = N_atoms.
n_bonds = 18

# N_rings = The cyclopentane ring.
n_rings = 1

# N_rotatable_bonds = The C-C bond connecting the ring to the carboxyl group,
# and the C-O single bond within the carboxyl group.
n_rotatable_bonds = 2

# Execute the calculation
calculate_bmc(product_name, n_atoms, n_bonds, n_rings, n_rotatable_bonds)

# Final answer for direct extraction
final_answer = (n_atoms * n_bonds * n_rings) / n_rotatable_bonds
# print(f'<<<{final_answer}>>>')