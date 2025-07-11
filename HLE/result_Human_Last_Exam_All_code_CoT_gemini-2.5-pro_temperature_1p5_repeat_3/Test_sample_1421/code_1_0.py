def display_fermionic_partition_function_formula():
    """
    This function prints the formula for the fermionic partition function Z
    in the imaginary time path integral formalism.
    """

    print("The formula for the fermionic partition function Z is given by a path integral over Grassmann fields:")
    print("-" * 80)

    # Main formula for Z
    print("\nZ = ∫ Dψ̄ Dψ exp(-S_E[ψ̄, ψ])\n")

    print("Where:")
    # Integral measure
    print("  - ∫ Dψ̄ Dψ: Represents the functional (or path) integral over all possible configurations")
    print("    of the independent Grassmann fields ψ̄(τ, x) and ψ(τ, x).")

    # Action S_E
    print("\n  - S_E[ψ̄, ψ]: Is the Euclidean action for the fermionic system, defined as:")
    print("    S_E = ∫_0^β dτ ∫ d^d x  L_E(ψ̄, ψ)")

    # Lagrangian L_E
    print("\n    Where L_E is the Euclidean Lagrangian density. For a system described by a")
    print("    single-particle Hamiltonian H₀, the Lagrangian is:")
    print("    L_E = ψ̄(τ, x) * [∂/∂τ + H₀] * ψ(τ, x)")

    # Components explained
    print("\n    In these expressions:")
    print("      * β = 1 / (k_B * T) is the inverse temperature.")
    print("      * τ is the imaginary time variable, integrated from 0 to β.")
    print("      * x represents the d spatial coordinates.")
    print("      * H₀ is the single-particle Hamiltonian operator (e.g., H₀ = -∇²/2m - μ).")
    print("      * ∂/∂τ is the partial derivative with respect to imaginary time.")

    # Boundary Conditions
    print("\n  - Boundary Conditions: The integral is performed subject to anti-periodic boundary")
    print("    conditions in the imaginary time direction, which is a fundamental property of fermions:")
    print("    ψ(τ = β, x) = -ψ(τ = 0, x)")
    print("    ψ̄(τ = β, x) = -ψ̄(τ = 0, x)")

    print("-" * 80)
    print("\nIn a more compact form, the partition function can also be expressed as the determinant of the operator in the Lagrangian:")
    print("\nZ = Det(∂/∂τ + H₀)")


# Execute the function to display the formula.
display_fermionic_partition_function_formula()
