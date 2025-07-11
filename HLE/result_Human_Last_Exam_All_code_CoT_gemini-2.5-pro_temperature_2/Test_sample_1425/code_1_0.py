def solve_partition_function():
    """
    This function presents the derivation and final formula for the partition function Z
    for the given system.
    """

    print("The partition function Z for a system with the grand canonical operator H = -μN is derived using the path integral for non-interacting bosons with zero single-particle energy.")
    print("The resulting equation is:")

    # Define the numerical components of the equation as per the instructions
    term_A = 1
    exponent_multiplier = -1

    # The final formula is constructed from these components and symbolic variables.
    # The full equation is Z = (term_A - exp(β*μ))**(exponent_multiplier * M)
    
    final_equation = "Z = (1 - exp(β*μ))**(-M)"
    print(f"\n{final_equation}\n")

    print("Where:")
    print("  Z: The Grand Canonical Partition Function")
    print("  M: The total number of available single-particle states")
    print("  β: The inverse temperature (1 / k_B * T)")
    print("  μ: The chemical potential")
    
    print("\nThe numbers present in the final equation are:")
    print(f"  - In the base (1 - exp(β*μ)), the number is: {term_A}")
    print(f"  - In the exponent (-M), the minus sign corresponds to the number: {exponent_multiplier}")

solve_partition_function()