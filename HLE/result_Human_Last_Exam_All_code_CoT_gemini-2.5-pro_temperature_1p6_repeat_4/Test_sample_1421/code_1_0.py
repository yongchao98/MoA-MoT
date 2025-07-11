def fermionic_partition_function_formula():
    """
    This function prints the formula for the fermionic partition function Z
    in the imaginary time representation using Feynman’s path integral formalism,
    along with an explanation of each term.
    """

    # --- Formula components ---
    z_equals = "Z ="
    integral_measure = "∫ Dψ̄ Dψ"
    exponential_part = "exp"
    action_integral = "- ∫₀ᵝ dτ ∫ dᵈr"
    lagrangian_density = "ψ̄(τ,r) [ ∂_τ - (ħ²/2m)∇² - μ ] ψ(τ,r)"

    # --- Print the full formula ---
    print("The formula for the fermionic partition function Z using the imaginary time path integral is:")
    # Using f-string for clear formatting of the equation
    print(f"\n  {z_equals} {integral_measure} {exponential_part} [ {action_integral} {lagrangian_density} ]\n")

    # --- Explanation of terms ---
    print("This formula is subject to anti-periodic boundary conditions: ψ(β,r) = -ψ(0,r)")
    print("\nWhere:")
    print("  Z: The Grand Canonical Partition Function.")
    print("  ∫ Dψ̄ Dψ: The functional path integral over all possible configurations of the fermionic fields.")
    print("  ψ(τ,r), ψ̄(τ,r): Anticommuting Grassmann fields, representing the fermions.")
    print("  τ: Imaginary time, which runs from 0 to β.")
    print("  β: Inverse temperature, defined as 1/(k_B T), where T is temperature and k_B is the Boltzmann constant.")
    print("  r: The spatial coordinate in d dimensions.")
    print("  ħ: The reduced Planck constant.")
    print("  m: The mass of the fermion.")
    print("  μ: The chemical potential.")
    print("  ∂_τ: The partial derivative with respect to imaginary time.")
    print("  ∇²: The Laplacian operator, representing the kinetic energy term.")
    print("\nThe entire expression in the exponent is the negative of the Euclidean action, S_E, for the system.")


if __name__ == '__main__':
    fermionic_partition_function_formula()
    # The final answer in the requested format
    final_formula = "Z = ∫ Dψ̄ Dψ exp[ - ∫₀ᵝ dτ ∫ dᵈr ψ̄(τ,r) ( ∂_τ - (ħ²/2m)∇² - μ ) ψ(τ,r) ] subject to anti-periodic boundary conditions ψ(β,r) = -ψ(0,r)"
    print(f"\n<<<{final_formula}>>>")