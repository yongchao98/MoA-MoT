import math

def display_partition_function():
    """
    This function prints the formula for the grand canonical partition function Z
    using the path integral formulation.
    """

    print("The grand canonical partition function Z for a system with Hamiltonian H' = H - μN, expressed using the path integral formalism, is as follows.")
    print("-" * 110)

    # Define the components of the formula as strings
    partition_function = "Z = ∫D[ψ*, ψ] * exp(-S[ψ*, ψ])"
    
    action_integral = "S[ψ*, ψ] = ∫[from 0 to β] dτ * ∫dᵈx * L_E"
    
    lagrangian_density = "L_E = ψ*(x, τ) * (∂/∂τ)ψ(x, τ) + H(ψ*, ψ) - μ*ψ*(x, τ)ψ(x, τ)"

    # Print the final equation with each part clearly defined
    print(f"The partition function is a functional integral over all field configurations:\n\n    {partition_function}\n")
    
    print("Where S[ψ*, ψ] is the Euclidean Action, which is the integral of the Euclidean Lagrangian density (L_E):")
    print(f"\n    {action_integral}\n")

    print("The Euclidean Lagrangian density L_E is given by:")
    print(f"\n    {lagrangian_density}\n")

    print("The symbols in the equation represent:")
    print(f"  - Z: The Grand Canonical Partition Function")
    print(f"  - ∫D[ψ*, ψ]: The functional ('path') integral over the fields ψ and its conjugate ψ*")
    print(f"  - S: The Euclidean Action")
    print(f"  - τ: Imaginary time, integrated from 0 to β")
    print(f"  - β: Inverse temperature, 1 / (k_B * T)")
    print(f"  - d: The number of spatial dimensions")
    print(f"  - H(ψ*, ψ): The original Hamiltonian density, with operators replaced by fields")
    print(f"  - μ: The chemical potential")
    print(f"  - N = ∫dᵈx ψ*ψ: The particle number")
    print("\nNote: The fields ψ have (anti-)periodic boundary conditions in τ: ψ(x, τ) = ±ψ(x, τ + β) for bosons/fermions.")
    print("-" * 110)

# Execute the function to display the result
display_partition_function()