def display_fermionic_partition_function():
    """
    This function prints the formula for the fermionic partition function Z
    in the imaginary time representation using Feynman's path integral formalism.
    The formula is built step-by-step with explanations.
    """

    print("The formula for the fermionic partition function Z is a path integral over Grassmann fields:")
    print("------------------------------------------------------------------------------------------\n")

    # Main formula for Z
    print("Z = ∫ Dψ̄ Dψ * exp(-S[ψ̄, ψ])\n")

    print("Where:")
    print("  - Z: The partition function.")
    print("  - ∫ Dψ̄ Dψ: A functional integral over all configurations of the anti-commuting")
    print("    Grassmann fields ψ(x, τ) and ψ̄(x, τ).")
    print("  - S[ψ̄, ψ]: The Euclidean action for the system.\n")

    # Formula for the action S
    print("The Euclidean action S is the integral of the Lagrangian density L over imaginary time τ and space:")
    print("---------------------------------------------------------------------------------------------------\n")
    print("  S[ψ̄, ψ] = ∫[from 0 to β] dτ ∫ dᵈx * L(ψ̄, ψ)\n")

    print("Where:")
    print("  - β: The inverse temperature, equal to 1 / (k_B * T). The integral runs over imaginary")
    print("    time τ from 0 to β.")
    print("  - ∫ dᵈx: An integral over d-dimensional space.")
    print("  - L: The Euclidean Lagrangian density. A common form is:\n")
    print("      L = ψ̄(x, τ) * [∂/∂τ + H₀] * ψ(x, τ)\n")
    print("      Here, H₀ is the single-particle Hamiltonian (e.g., H₀ = -∇²/ (2 * m) - μ), where μ is the chemical potential.\n")


    # Boundary Conditions
    print("Crucially, for fermions, the integral is over fields that satisfy anti-periodic boundary conditions:")
    print("-----------------------------------------------------------------------------------------------------\n")
    print("  ψ(x, τ = β) = -ψ(x, τ = 0)")
    print("  ψ̄(x, τ = β) = -ψ̄(x, τ = 0)\n")

# Execute the function to display the formula.
display_fermionic_partition_function()