def display_fermionic_partition_function_formula():
    """
    This script prints the formula for the fermionic partition function (Z)
    in the imaginary time path integral formalism.
    """

    print("The formula for the fermionic partition function Z in the imaginary time path integral formalism is:")
    print("\n" + "="*80 + "\n")

    # The main formula is constructed and printed piece by piece.
    print("Z = ∫ D[ψ*(τ), ψ(τ)] * exp(-S[ψ*, ψ])")

    print("\n" + "="*80 + "\n")
    print("Where the components of the formula are defined as follows:\n")

    # 1. Z: The Partition Function
    print("Z: The partition function of the fermionic system in thermal equilibrium.")

    # 2. ∫ D[ψ*, ψ(τ)]: The Functional Integral
    print("\n∫ D[ψ*(τ), ψ(τ)]: The Feynman path integral.")
    print("  - This integrates over all possible configurations (paths) of the fermionic fields.")
    print("  - ψ(τ) and ψ*(τ) are independent Grassmann fields. These are anti-commuting numbers")
    print("    that mathematically represent the Pauli exclusion principle for fermions.")
    print("  - τ is the imaginary time variable, which runs on the interval [0, β].")

    # 3. S[ψ*, ψ]: The Euclidean Action
    print("\nS[ψ*, ψ]: The Euclidean action for the fermionic fields.")
    print("S = ∫_0^β dτ [ ψ*(τ) * (∂/∂τ)ψ(τ) + H(ψ*(τ), ψ(τ)) ]")
    print("  - ∫_0^β dτ: Integration over the imaginary time interval, where β = 1 / (k_B * T)")
    print("    (T is temperature, k_B is the Boltzmann constant).")
    print("  - H(ψ*(τ), ψ(τ)): The system's Hamiltonian, where creation and annihilation operators")
    print("    have been replaced by the corresponding Grassmann fields ψ* and ψ.")

    # 4. Anti-Periodic Boundary Conditions
    print("\nBoundary Conditions: The path integral for fermions requires ANTI-PERIODIC boundary conditions")
    print("in imaginary time. This is a crucial feature that distinguishes them from bosons.")
    print("The fields must satisfy:")
    print("  ψ(τ + β) = -ψ(τ)")
    print("  ψ*(τ + β) = -ψ*(τ)")

    print("\n" + "="*80)


if __name__ == "__main__":
    display_fermionic_partition_function_formula()
