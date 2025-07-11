def print_fermionic_partition_function_formula():
    """
    Prints the formula for the fermionic partition function Z in the
    imaginary time path integral formalism.
    """

    # Define the components of the formula using unicode for mathematical symbols
    partition_function = "Z"
    integral_sign = "∫"
    differential_measure = "Dψ̄ Dψ"
    exponential_part = "exp(-S_E[ψ̄, ψ])"

    # Assemble the main formula
    formula = f"{partition_function} = {integral_sign} {differential_measure} {exponential_part}"

    # Define the Euclidean action
    action = "S_E[ψ̄, ψ]"
    action_definition = "∫₀ᵝ dτ ∫ dᵈx L_E(ψ̄, ψ, ∂ψ̄, ∂ψ)"

    # Define the boundary conditions
    boundary_condition = "ψ(τ, x) = -ψ(τ + β, x)"
    conjugate_boundary_condition = "ψ̄(τ, x) = -ψ̄(τ + β, x)"

    # Print the full explanation
    print("The formula for the fermionic partition function Z is:")
    print("-" * 50)
    print(formula)
    print("-" * 50)
    print("\nWhere:")
    print(f"  - Z: The partition function.")
    print(f"  - {integral_sign} {differential_measure}: The functional integral (Berezin integral) over all anti-commuting Grassmann field configurations ψ and ψ̄.")
    print(f"  - S_E: The Euclidean action for the fermionic fields, defined as:")
    print(f"      {action} = {action_definition}")
    print(f"  - L_E: The Euclidean Lagrangian density.")
    print(f"  - τ: Imaginary time, integrated from 0 to β.")
    print(f"  - β: Inverse temperature, β = 1/(k_B T).")
    print(f"\nThe integral is performed subject to anti-periodic boundary conditions in imaginary time:")
    print(f"  - {boundary_condition}")
    print(f"  - {conjugate_boundary_condition}")

if __name__ == "__main__":
    print_fermionic_partition_function_formula()