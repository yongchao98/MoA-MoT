def calculate_molecular_weight():
    """
    Calculates the molecular weight of 4-pentyl-benzonitrile (C12H15N)
    based on the provided design specifications.
    """
    # Standard atomic weights (g/mol)
    atomic_weight_C = 12.011
    atomic_weight_H = 1.008
    atomic_weight_N = 14.007

    # The molecule is CnH2n+1-Ph-CN, where Ph is a benzene ring (C6H4)
    # From the prompt, we start with a pentyl chain, so n=5.
    n = 5

    # Calculate the number of atoms for each element
    # Carbons: n from alkyl chain + 6 from benzene ring + 1 from cyano group
    num_C = n + 6 + 1
    # Hydrogens: (2n+1) from alkyl chain + 4 from benzene ring
    num_H = (2 * n + 1) + 4
    # Nitrogens: 1 from the cyano group
    num_N = 1

    # Calculate the contribution of each element to the total molecular weight
    mass_C = num_C * atomic_weight_C
    mass_H = num_H * atomic_weight_H
    mass_N = num_N * atomic_weight_N

    # Calculate the total molecular weight
    total_mw = mass_C + mass_H + mass_N

    # Print the details of the calculation as an equation
    print("Calculating the molecular weight for 4-pentyl-benzonitrile (C12H15N):")
    print(f"Formula: C{num_C}H{num_H}N{num_N}")
    print("\nEquation:")
    # The instruction is to output each number in the final equation.
    print(f"({num_C} * {atomic_weight_C}) + ({num_H} * {atomic_weight_H}) + ({num_N} * {atomic_weight_N}) = {total_mw:.3f} g/mol")

    # Print the final result clearly
    print(f"\nThe calculated molecular weight is: {total_mw:.3f} g/mol")

# Execute the function
calculate_molecular_weight()