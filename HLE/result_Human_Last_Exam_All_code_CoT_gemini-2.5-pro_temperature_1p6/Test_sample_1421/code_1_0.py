def print_fermionic_partition_function_formula():
    """
    This function prints the formula for the fermionic partition function (Z)
    in the imaginary time representation using Feynman's path integral formalism.
    """

    # --- Main Formula for the Partition Function Z ---
    partition_function = "Z = ∫ D[ψ_bar] D[ψ] * exp(-S_E)"
    
    # --- Euclidean Action S_E ---
    # The action S_E is defined by an integral from 0 to β
    euclidean_action = "S_E = ∫_0^β dτ ∫ d^d x * L_E(ψ_bar, ψ)"
    
    # --- Boundary Conditions ---
    # The boundary conditions are anti-periodic between τ=0 and τ=β
    boundary_condition = "ψ(x, 0) = -ψ(x, β)"

    # --- Print the complete formula with explanations ---
    print("The formula for the fermionic partition function Z is given by:")
    print("----------------------------------------------------------------")
    
    # Print Partition Function
    print("\n1. Partition Function (Z):")
    print(f"   {partition_function}")
    print("\n   Where:")
    print("   - ∫ D[ψ_bar] D[ψ] is the path integral over all configurations of the fermionic Grassmann fields ψ and ψ_bar.")
    print("   - S_E is the Euclidean action of the system.")
    
    # Print Euclidean Action and highlight the numbers
    print("\n2. Euclidean Action (S_E):")
    print(f"   {euclidean_action}")
    print("\n   Where:")
    print("   - L_E is the Euclidean Lagrangian density.")
    print("   - The first integral is over imaginary time (τ) from the number 0 to β.")
    print("   - The second integral is over the d spatial dimensions (x).")
    
    # Print Boundary Conditions and highlight the numbers
    print("\n3. Boundary Conditions:")
    print("   The path integral is subject to anti-periodic boundary conditions in imaginary time:")
    print(f"   {boundary_condition}")
    print("\n   Where:")
    print(f"   - The condition relates the field at imaginary time τ=0 to its value at τ=β.")
    
    # Print Definition of Beta
    print("\n4. Definition of β:")
    print("   β = 1 / (k_B * T), where T is temperature and k_B is the Boltzmann constant.")
    print("----------------------------------------------------------------")


# Execute the function to display the formula
print_fermionic_partition_function_formula()