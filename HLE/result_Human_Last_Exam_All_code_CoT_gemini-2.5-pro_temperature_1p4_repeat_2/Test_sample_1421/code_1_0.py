def print_fermionic_partition_function_formula():
    """
    Prints the formula for the fermionic partition function Z in the
    imaginary time representation using Feynman's path integral formalism.
    """

    # Define the components of the formula
    Z = "Z"
    integral_sign = "∫"
    measure = "D[ψ̄] D[ψ]"
    exponent = "exp(-S[ψ̄, ψ])"
    # Construct the final equation string
    final_equation = f"{Z} = {integral_sign} {measure} {exponent}"

    action_S = "S[ψ̄, ψ]"
    integral_time = "∫ from 0 to β dτ"
    lagrangian = "L(ψ̄(τ), ψ(τ))"
    action_definition = f"{action_S} = {integral_time} {lagrangian}"

    lagrangian_example = "L = ψ̄(τ) (∂_τ + H) ψ(τ)"

    boundary_condition = "ψ(β) = -ψ(0)"

    # Print the explanation and the formulas
    print("The formula for the fermionic partition function Z in the imaginary time path integral formalism is as follows.")
    print("="*80)
    print("Final Equation:")
    print(f"{Z} = {integral_sign} {measure} {exponent}")
    print("="*80)

    print("\nWhere the terms are defined as:\n")
    print(f"*   Z: The grand canonical partition function of the fermionic system.")
    print(f"*   ∫ D[ψ̄] D[ψ]: The functional integral (or Berezin integral) over all possible configurations of the Grassmann fields ψ̄ and ψ.")
    print(f"*   ψ(τ) and ψ̄(τ): Independent Grassmann-valued fields, which are functions of imaginary time τ. They represent the fermionic degrees of freedom.")
    print(f"*   The fields must satisfy the anti-periodic boundary condition in imaginary time: {boundary_condition}, where β = 1/(k_B T) is the inverse temperature.")
    print(f"*   S[ψ̄, ψ]: The Euclidean action of the system, which is defined as the integral of the Euclidean Lagrangian (L) over imaginary time:")
    print(f"\n      {action_definition}\n")
    print(f"*   L: The Euclidean Lagrangian. For a system described by a single-particle Hamiltonian H (which includes kinetic energy and chemical potential), a common form is:")
    print(f"\n      {lagrangian_example}\n")

if __name__ == '__main__':
    print_fermionic_partition_function_formula()
