def print_fermionic_partition_function_formula():
    """
    This function prints the formula for the fermionic partition function Z
    in the imaginary time path integral formalism, along with explanations
    for each component.
    """

    # --- Main Formulas ---
    formula_header = "The formula for the fermionic partition function (Z) using the path integral in imaginary time is:"
    formula_Z = "Z = ∫ D[ψ*, ψ] exp(-S[ψ*, ψ])"

    action_header = "\nThe action (S) is the integral of the Lagrangian (L) over imaginary time from 0 to β:"
    formula_S = "S[ψ*, ψ] = ∫[from τ=0 to β] dτ L(ψ*, ψ, τ)"

    lagrangian_header = "\nA common form for the Lagrangian (L) is:"
    formula_L = "L(ψ*, ψ, τ) = ψ*(τ)∂_τψ(τ) + H(ψ*, ψ)"

    boundary_condition_header = "\nThe integral is performed subject to anti-periodic boundary conditions in imaginary time:"
    boundary_condition = "ψ(β) = -ψ(0)   and   ψ*(β) = -ψ*(0)"

    # --- Explanation of Terms ---
    explanation_header = "\n" + "="*60 + "\nExplanation of Terms:\n" + "="*60
    explanation_Z = "Z: The fermionic partition function."
    explanation_integral = "∫ D[ψ*, ψ]: The functional or path integral over all possible field configurations."
    explanation_psi = "ψ(τ), ψ*(τ): Anti-commuting Grassmann fields, which represent the fermions. They are functions of imaginary time τ."
    explanation_S = "S[ψ*, ψ]: The Euclidean action of the system."
    explanation_beta = "β: Inverse temperature, defined as β = 1 / (k_B * T), where k_B is the Boltzmann constant and T is the temperature."
    explanation_tau = "τ: Imaginary time, which runs from 0 to β."
    explanation_L = "L(ψ*, ψ, τ): The Lagrangian of the system in imaginary time."
    explanation_d_tau = "∂_τ: The partial derivative with respect to imaginary time."
    explanation_H = "H(ψ*, ψ): The Hamiltonian of the system, which describes its energy (e.g., kinetic and potential/interaction terms)."
    explanation_bc = "Anti-periodic boundary conditions: This is a key feature for fermions, arising from the trace operation in the definition of Z."
    explanation_numbers = "The numbers in the final equation are the integration limits for the action: 0 and β."

    # --- Print Everything ---
    print(formula_header)
    print(f"  {formula_Z}")

    print(action_header)
    print(f"  {formula_S}")

    print(lagrangian_header)
    print(f"  {formula_L}")

    print(boundary_condition_header)
    print(f"  {boundary_condition}")

    print(explanation_header)
    print(f"- {explanation_Z}")
    print(f"- {explanation_integral}")
    print(f"- {explanation_psi}")
    print(f"- {explanation_S}")
    print(f"- {explanation_beta}")
    print(f"- {explanation_tau}")
    print(f"- {explanation_L}")
    print(f"- {explanation_d_tau}")
    print(f"- {explanation_H}")
    print(f"- {explanation_bc}")
    print(f"- {explanation_numbers}")

# Execute the function to print the formula and explanation.
if __name__ == "__main__":
    print_fermionic_partition_function_formula()