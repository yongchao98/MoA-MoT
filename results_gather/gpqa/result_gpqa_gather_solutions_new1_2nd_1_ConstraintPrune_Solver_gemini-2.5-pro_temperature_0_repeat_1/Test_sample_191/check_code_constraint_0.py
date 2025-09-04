import sympy

def check_answer():
    """
    Checks the correctness of the LLM's answer based on physics principles.
    
    The core principle is electrostatic shielding: the external electric field
    cannot depend on the internal geometry of the cavity (l, s, theta).
    It should only depend on the distance from the center of the conductor (L).
    """
    
    # Define symbolic variables
    q, epsilon_o, pi = sympy.symbols('q epsilon_o pi')
    L, l, s, theta = sympy.symbols('L l s theta')

    # The options provided in the question
    # Note: The constant factor 1/(4*pi*epsilon_o) is common to all and can be ignored for dependency checking.
    options = {
        'A': q / L**2,
        'B': q / (l - s * sympy.cos(theta))**2,
        'C': q / l**2,
        'D': q / (l + s * sympy.cos(theta))**2
    }

    # The final answer provided by the LLM
    llm_answer = 'A'

    # --- Verification Logic ---

    # 1. Check if the provided answer is one of the options
    if llm_answer not in options:
        return f"Incorrect. The answer '{llm_answer}' is not one of the valid options A, B, C, or D."

    # 2. Get the formula for the chosen answer
    chosen_formula = options[llm_answer]
    
    # 3. Define the set of variables related to the internal geometry (forbidden variables)
    forbidden_vars = {l, s, theta}
    
    # 4. Check if the chosen formula violates the shielding principle
    if any(var in chosen_formula.free_symbols for var in forbidden_vars):
        violating_vars = {str(var) for var in forbidden_vars if var in chosen_formula.free_symbols}
        return (f"Incorrect. The answer {llm_answer} corresponds to the formula E ∝ {chosen_formula}, "
                f"which is wrong because it violates the principle of electrostatic shielding. "
                f"The external electric field cannot depend on the internal geometry variables: {', '.join(violating_vars)}.")

    # 5. Check if the chosen formula depends on the correct variable L
    if L not in chosen_formula.free_symbols:
        return (f"Incorrect. The answer {llm_answer} corresponds to the formula E ∝ {chosen_formula}, "
                f"which is wrong. The external electric field must depend on the distance 'L' from the "
                f"center of the spherical conductor.")

    # 6. Verify that all other options are incorrect
    for option_key, formula in options.items():
        if option_key == llm_answer:
            continue
        # An incorrect option should violate the shielding principle
        if not any(var in formula.free_symbols for var in forbidden_vars):
            return (f"Check failed. The logic suggests that option {option_key} might also be correct, "
                    f"as it does not depend on the forbidden variables {forbidden_vars}. The check is inconclusive.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_answer()
print(result)