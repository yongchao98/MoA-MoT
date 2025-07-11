def display_partition_function():
    """
    This function prints the path integral representation of the grand
    canonical partition function Z for a Hamiltonian H' = H - μN.
    """

    # Define the components of the final equation as strings
    partition_function_symbol = "Z"
    equals = "="
    integral_measure = "∫ D[ψ*, ψ]"
    exponential_term = "exp"
    
    # Define the Action S
    action_symbol = "S"
    action_definition = "∫₀^β dτ ∫ dᵈx [ ψ*(x,τ)ħ∂_τψ(x,τ) + H(ψ*,ψ) - μN(ψ*,ψ) ]"
    
    # Header
    print("The grand canonical partition function Z is given by the path integral:")
    print("-" * 70)

    # Print the general form
    print("General form: Z = ∫ D[ψ*, ψ] * exp(-S)")
    print("Where S is the Euclidean action.\n")

    # Print the detailed final equation
    print("The full expression for the partition function Z is:")
    
    # This line constructs and prints the final equation, showing each symbolic component.
    final_equation = (
        f"{partition_function_symbol} {equals} {integral_measure} "
        f"{exponential_term}{{ -{action_definition} }}"
    )
    print(final_equation)
    print("-" * 70)
    
    # Explanation of terms
    print("Where:")
    print("  Z      : The grand canonical partition function.")
    print("  ∫D[ψ*,ψ]: The functional (path) integral over all field configurations.")
    print("  ψ, ψ*  : The fields (Grassmann variables for fermions, complex for bosons).")
    print("  β      : Inverse temperature (1 / k_B T).")
    print("  ħ      : The reduced Planck constant.")
    print("  τ      : Imaginary time, integrated from 0 to β.")
    print("  x, dᵈx : Position vector and volume element in d spatial dimensions.")
    print("  ∂_τ    : The partial derivative with respect to imaginary time.")
    print("  H(ψ*,ψ): The Hamiltonian written as a function of the fields.")
    print("  μ      : The chemical potential.")
    print("  N(ψ*,ψ): The particle number operator written as a function of the fields.")
    print("\nThe fields ψ must satisfy periodic (for bosons) or anti-periodic (for fermions) boundary conditions in imaginary time: ψ(x, 0) = ±ψ(x, β).")

if __name__ == "__main__":
    display_partition_function()