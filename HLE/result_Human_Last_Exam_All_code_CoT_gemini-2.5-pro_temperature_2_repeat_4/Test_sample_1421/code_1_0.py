def print_fermionic_partition_function_formula():
    """
    This function prints the formula for the fermionic partition function Z
    in the imaginary time path integral formalism.
    """

    # Title
    print("Formula for the Fermionic Partition Function (Z)")
    print("=" * 50)

    # Main formula for Z
    partition_function = "Z = ∫ D[ψ*, ψ] exp(-S[ψ*, ψ])"
    print("\nThe partition function Z is given by the functional integral:")
    print(f"  {partition_function}")

    # Explanation of Z
    print("\nWhere:")
    print("  - '∫ D[ψ*, ψ]' is the functional integral over all configurations of the")
    print("    Grassmann fields ψ(τ) and ψ*(τ).")
    print("  - The integral is taken over fields that satisfy anti-periodic")
    print("    boundary conditions in imaginary time: ψ(β) = -ψ(0).")
    print("  - S[ψ*, ψ] is the Euclidean action.")

    # Formula for the Euclidean Action S
    action = "S = ∫₀ᵝ dτ [ Σₖ ψₖ*(τ) ∂_τ ψₖ(τ) + H(ψ*(τ), ψ(τ)) ]"
    print("\nThe Euclidean action S for the fermions is defined as:")
    print(f"  {action}")

    # Explanation of S
    print("\nWhere:")
    print("  - '∫₀ᵝ dτ' is the integral over imaginary time τ from 0 to β.")
    print("  - 'β' is the inverse temperature, defined as β = 1 / (k_B * T).")
    print("  - 'ψₖ(τ)' and 'ψₖ*(τ)' are anticommuting Grassmann fields, which represent")
    print("    fermions in the path integral. The index 'k' labels the")
    print("    single-particle states (e.g., momentum and spin).")
    print("  - '∂_τ' represents the derivative with respect to imaginary time.")
    print("  - 'H(ψ*, ψ)' is the system's Hamiltonian, where the original fermion")
    print("    operators have been replaced by the corresponding Grassmann fields.")

if __name__ == "__main__":
    print_fermionic_partition_function_formula()