def print_fermionic_partition_function_formula():
    """
    This function prints the formula for the fermionic partition function Z
    in the imaginary time path integral formalism.
    """

    # Introduction to the formula
    print("The formula for the fermionic partition function Z is given by a functional integral over anti-commuting Grassmann fields.")
    print("-" * 80)

    # 1. The main formula for Z
    print("1. The Partition Function (Z):")
    # Using psi_bar for the conjugate field psi-bar
    # Using D[psi_bar]D[psi] for the functional integration measure
    # Using S_E for the Euclidean action
    print("   Z = ∫ D[ψ_bar] D[ψ] * exp(-S_E[ψ_bar, ψ])")
    print("\n   Where:")
    print("   - ∫ D[ψ_bar] D[ψ] is the functional integral over all possible field configurations.")
    print("   - ψ(τ, x) and ψ_bar(τ, x) are the Grassmann-valued fermion fields.")
    print("   - S_E is the Euclidean action of the system.")
    print("-" * 80)

    # 2. The Euclidean Action (S_E)
    print("2. The Euclidean Action (S_E):")
    # Using beta for the inverse temperature 1/(k_B * T)
    # Using tau for imaginary time
    # Using L_E for the Euclidean Lagrangian
    print("   S_E = ∫[from 0 to β] dτ L_E(τ)")
    print("\n   Where:")
    print("   - The integral is over imaginary time 'τ' from 0 to β.")
    print("   - β = 1 / (k_B * T) is the inverse temperature.")
    print("   - L_E is the Euclidean Lagrangian.")
    print("-" * 80)

    # 3. The Euclidean Lagrangian (L_E)
    print("3. The Euclidean Lagrangian (L_E):")
    # Using H_0 for the single-particle Hamiltonian
    # Using ∂_τ for the partial derivative with respect to tau
    print("   L_E = ∫ d³x ψ_bar(τ, x) * [ħ∂_τ + H_0] * ψ(τ, x)")
    print("\n   Where:")
    print("   - The integral is over all spatial dimensions 'x'.")
    print("   - ħ is the reduced Planck constant.")
    print("   - ∂_τ is the partial derivative with respect to imaginary time.")
    print("   - H_0 is the single-particle Hamiltonian (e.g., H_0 = -ħ²∇²/2m - μ).")
    print("-" * 80)

    # 4. The Boundary Conditions
    print("4. Anti-Periodic Boundary Conditions:")
    print("   A key feature for fermions is that the Grassmann fields must satisfy anti-periodic boundary conditions in imaginary time:")
    print("   ψ(τ + β, x) = -ψ(τ, x)")
    print("   ψ_bar(τ + β, x) = -ψ_bar(τ, x)")
    print("\n   This arises from the anti-commutation relations of fermion operators when taking the trace.")
    print("-" * 80)

if __name__ == '__main__':
    print_fermionic_partition_function_formula()