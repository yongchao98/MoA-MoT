def display_fermionic_partition_function():
    """
    This function prints the formula for the fermionic partition function Z
    in the imaginary time path integral formalism.
    """

    # Define the components of the formula
    z_definition = "Z = ∫ D[ψ̄] D[ψ] * exp(-S_E[ψ̄, ψ])"
    se_definition = "S_E[ψ̄, ψ] = ∫_0^β dτ ∫ d^d x  L_E(ψ̄(τ, x), ψ(τ, x))"
    bc_definition = "ψ(τ, x) = -ψ(τ + β, x)  and  ψ̄(τ, x) = -ψ̄(τ + β, x)"
    le_example = "L_E = ψ̄(τ, x) * [ħ ∂/∂τ - (ħ²/2m)∇² - μ] * ψ(τ, x)"

    # Print the explanation and the formulas
    print("The formula for the fermionic partition function (Z) in the imaginary time representation is:")
    print("-" * 80)
    print(f"Z = ∫ D[ψ̄] D[ψ] * exp( - ∫_0^β dτ ∫ d^d x  L_E )")
    print("-" * 80)
    print("This formula is composed of the following key parts:\n")

    print("1. The Path Integral:")
    print("   The integral ∫ D[ψ̄] D[ψ] is a functional integral over all configurations of")
    print("   the anti-commuting Grassmann fields ψ and ψ̄.\n")

    print("2. The Exponential Weight:")
    print("   The term exp(-S_E) weights each path by the Euclidean action S_E.\n")

    print("3. The Euclidean Action (S_E):")
    print("   This is the integral of the Euclidean Lagrangian density (L_E) over imaginary")
    print(f"   time τ from 0 to β = 1/(k_B*T) and all spatial dimensions (d^d x).")
    print(f"   S_E = ∫_0^β dτ ∫ d^d x  L_E\n")
    print("   For a non-relativistic fermion, a common example for L_E is:")
    print(f"   {le_example}\n")

    print("4. The Boundary Conditions:")
    print("   The integral is restricted to fields that are anti-periodic in imaginary time:")
    print(f"   {bc_definition}\n")

if __name__ == '__main__':
    display_fermionic_partition_function()