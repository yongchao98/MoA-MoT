def print_fermionic_partition_function():
    """
    Prints the formula for the fermionic partition function Z in the
    imaginary time path integral formalism and explains its components.
    """

    # The formula is presented using common physics notation.
    # Z = The partition function
    # ∫ D[ψ_bar, ψ] = The functional integral over Grassmann fields
    # e^(-S) = The Boltzmann factor with the Euclidean action S
    formula_Z = "Z = ∫ D[ψ_bar, ψ] * e^(-S[ψ_bar, ψ])"

    # S = The Euclidean action
    # ∫dτ from 0 to β = Integral over imaginary time
    # ∫d^d x = Integral over spatial dimensions
    # L_E = Euclidean Lagrangian density
    formula_S = "S[ψ_bar, ψ] = ∫_0^β dτ ∫d^d x L_E(ψ_bar, ψ)"

    # L_E = Euclidean Lagrangian density
    # ψ_bar, ψ = Fermionic Grassmann fields
    # ∂_τ = Partial derivative with respect to imaginary time
    # H_0 = Single-particle Hamiltonian (e.g., kinetic energy + chemical potential)
    formula_L = "L_E = ψ_bar * (∂_τ + H_0) * ψ"
    formula_L_example = "Example H_0: H_0 = - (ħ²/2m)∇² - μ"


    # Anti-periodic boundary conditions for fermions
    boundary_conditions = "Boundary Conditions: ψ(β, x) = -ψ(0, x)"

    print("The formula for the fermionic partition function (Z) is:")
    print("=" * 50)
    print(formula_Z)
    print("-" * 50)

    print("\nWhere the components of the formula are:\n")

    print("1. The Euclidean Action (S):")
    print(f"   {formula_S}")
    print("\n   With a typical Euclidean Lagrangian (L_E) for non-relativistic fermions:")
    print(f"   {formula_L}")
    print(f"   where a common single-particle Hamiltonian (H_0) is: {formula_L_example}")


    print("\n2. Boundary Conditions:")
    print("   The path integral is calculated with anti-periodic boundary conditions in imaginary time, which is a defining feature of fermions.")
    print(f"   {boundary_conditions}")


    print("\n3. Glossary of Symbols:")
    print("   Z: The partition function.")
    print("   ∫ D[ψ_bar, ψ]: The functional (or path) integral over all possible configurations of the fermionic fields.")
    print("   ψ(τ, x), ψ_bar(τ, x): Anti-commuting Grassmann fields representing the fermions.")
    print("   S[ψ_bar, ψ]: The Euclidean action in imaginary time.")
    print("   τ: Imaginary time, integrated from 0 to β.")
    print("   β: Inverse temperature, defined as 1 / (k_B * T).")
    print("   k_B: Boltzmann constant.")
    print("   T: Temperature.")
    print("   x: Spatial coordinate(s), integrated over the volume of the system.")
    print("   d: Number of spatial dimensions.")
    print("   ħ: Reduced Planck constant.")
    print("   m: Mass of the fermion.")
    print("   μ: Chemical potential.")
    print("   ∇²: The Laplacian operator (spatial second derivative).")
    print("=" * 50)

# Execute the function to print the formula
print_fermionic_partition_function()