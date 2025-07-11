def print_fermionic_partition_function_formula():
    """
    Prints the formula for the fermionic partition function Z in the
    imaginary time path integral formalism.
    """

    # Define the components of the formula as strings
    partition_function = "Z"
    equals_sign = "="
    integral_symbol = "∫"
    path_integral_measure = "D[ψ*] D[ψ]"
    exponential_term = "e^(-S[ψ*, ψ])"
    
    action = "S[ψ*, ψ]"
    action_definition = "∫_0^β dτ ∫ d^d x  ψ*(τ, x) [∂_τ + H] ψ(τ, x)"
    
    boundary_condition = "ψ(β, x) = -ψ(0, x)"

    # Print the main formula
    print("The formula for the fermionic partition function is:")
    print(f"  {partition_function} {equals_sign} {integral_symbol} {path_integral_measure} {exponential_term}\n")
    
    # Print the definitions of the terms
    print("Where:")
    print(f"- {partition_function}: The Partition Function.")
    print(f"- {integral_symbol} {path_integral_measure}: The path integral over all anti-commuting Grassmann field configurations ψ(τ, x) and ψ*(τ, x).")
    print(f"- {exponential_term}: The Boltzmann weight, with S being the Euclidean action.\n")

    print("The Euclidean action S is given by:")
    print(f"  {action} {equals_sign} {action_definition}\n")

    print("And the components of the action are:")
    print("- τ: Imaginary time, integrated from 0 to β.")
    print("- β: Inverse temperature (1 / k_B T).")
    print("- x: Spatial coordinates.")
    print("- ψ(τ, x), ψ*(τ, x): The Grassmann fields representing the fermions.")
    print("- ∂_τ: The partial derivative with respect to imaginary time.")
    print("- H: The single-particle Hamiltonian operator (e.g., kinetic energy + potential).\n")

    print("A crucial element for fermions is the anti-periodic boundary condition in imaginary time:")
    print(f"  {boundary_condition}\n")

if __name__ == '__main__':
    print_fermionic_partition_function_formula()
