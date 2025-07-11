def calculate_bottcher_complexity():
    """
    Calculates the Böttcher Molecular Complexity for the product of the
    Favorskii rearrangement of 2-chlorocyclohexanone.
    """
    # The product is cyclopentanecarboxylic acid (C6H10O2).
    
    # N_A: Number of non-hydrogen atoms
    num_carbon = 6
    num_oxygen = 2
    N_A = num_carbon + num_oxygen

    # N_T: Total number of atoms
    num_hydrogen = 10
    N_T = num_carbon + num_oxygen + num_hydrogen

    # N_B: Number of bonds between non-hydrogen atoms (sum of bond orders).
    # The structure has:
    # - 5 C-C single bonds in the cyclopentane ring
    # - 1 C-C single bond connecting the ring to the carboxyl group
    # - 1 C=O double bond
    # - 1 C-O single bond
    N_B = 5 * 1 + 1 * 1 + 1 * 2 + 1 * 1

    # Calculate the Böttcher Molecular Complexity (BMC)
    # BMC = (N_A * N_B) / (N_T - 1)
    if (N_T - 1) == 0:
        bmc_value = float('inf')
    else:
        bmc_value = (N_A * N_B) / (N_T - 1)

    # Print the explanation, the final equation with its components, and the result.
    print("The product of the reaction is cyclopentanecarboxylic acid (C6H10O2).")
    print("The Böttcher Molecular Complexity (BMC) is calculated using the formula:")
    print("BMC = (N_A * N_B) / (N_T - 1)\n")
    print("Where for cyclopentanecarboxylic acid:")
    print(f"  N_A (non-hydrogen atoms) = {N_A}")
    print(f"  N_B (bonds between non-hydrogen atoms) = {N_B}")
    print(f"  N_T (total atoms) = {N_T}\n")
    
    print("Final Calculation:")
    # The prompt requires printing each number in the final equation
    print(f"BMC = ({N_A} * {N_B}) / ({N_T} - 1)")
    print(f"BMC = {N_A * N_B} / {N_T - 1}")
    print(f"BMC = {bmc_value}")

if __name__ == "__main__":
    calculate_bottcher_complexity()