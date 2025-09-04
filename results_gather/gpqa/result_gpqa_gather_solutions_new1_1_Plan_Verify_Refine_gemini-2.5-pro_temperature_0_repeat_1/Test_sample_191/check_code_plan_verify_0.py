import sympy

def check_answer_correctness():
    """
    Checks the correctness of the answer based on the principle of electrostatic shielding.

    The electric field outside a spherical conductor with a cavity depends only on the
    total charge induced on the outer surface and the distance from the conductor's center (L).
    It must be independent of the internal geometry (s, l, theta).
    """
    # Define all variables symbolically to analyze the formulas
    # We only need to know what the variables are, not their values.
    q, L, l, s, theta, epsilon_o = sympy.symbols('q L l s theta epsilon_o')

    # The constant k = 1 / (4 * pi * epsilon_o) is common to all options and doesn't affect dependencies.
    k = 1

    # Define the formulas for each option
    options = {
        'A': k * q / (l - s * sympy.cos(theta))**2,
        'B': k * q / (l + s * sympy.cos(theta))**2,
        'C': k * q / L**2,
        'D': k * q / l**2
    }

    # The final answer provided by the LLM to be checked
    llm_answer = 'C'

    # --- Verification Logic ---

    # Principle: The correct formula must not depend on internal geometry.
    # These variables describe the internal geometry and should NOT be in the final formula.
    internal_geometry_vars = {s, l, theta}

    # Principle: The formula must depend on the distance from the conductor's center.
    external_distance_var = L

    # Find the option that satisfies the physical principles
    valid_option = None
    reasons_for_rejection = {}

    for option_key, formula in options.items():
        # Get all variables present in the formula
        variables_in_formula = formula.free_symbols

        # Check for violation of electrostatic shielding
        violating_vars = internal_geometry_vars.intersection(variables_in_formula)
        if violating_vars:
            reasons_for_rejection[option_key] = (
                f"violates the principle of electrostatic shielding by incorrectly "
                f"depending on internal geometry parameter(s): {violating_vars}."
            )
            continue

        # Check if it correctly depends on the external distance L
        if external_distance_var not in variables_in_formula:
            reasons_for_rejection[option_key] = (
                f"incorrectly omits the dependency on the distance 'L' from the "
                f"conductor's center."
            )
            continue

        # If no violations, this is the physically correct option
        valid_option = option_key

    # --- Final Verdict ---
    if valid_option is None:
        return "Error in checking logic or problem statement: No valid option found."

    if llm_answer == valid_option:
        return "Correct"
    else:
        reason = reasons_for_rejection.get(llm_answer, "it is not the correct formula derived from physical principles.")
        return (f"Incorrect. The provided answer is '{llm_answer}', which is wrong because it {reason}\n"
                f"The correct answer is '{valid_option}', as it is the only option that correctly "
                f"depends only on the total charge 'q' and the distance 'L' from the conductor's center, "
                f"adhering to the principle of electrostatic shielding.")

# Run the check
result = check_answer_correctness()
print(result)