import sympy

def check_larmor_precession_answer():
    """
    Checks the correctness of the answer for the Larmor precession problem.

    The problem asks for the oscillation frequency of the z-component of a spin-1/2
    particle's magnetic moment when a magnetic field is suddenly switched from the
    Z-direction to the Y-direction.

    The correct physical principle is Larmor precession, and the oscillation frequency
    (Larmor frequency) is given by ω = γB.
    """

    # Define symbolic variables to represent the physical quantities
    gamma, B = sympy.symbols('gamma B', real=True, positive=True)

    # Define the options provided in the question as symbolic expressions
    options = {
        'A': gamma * B / 2,
        'B': sympy.sqrt(2) * gamma * B,
        'C': gamma * B,
        'D': gamma * B / sympy.sqrt(2)
    }

    # The correct answer based on the physics of Larmor precession
    correct_physical_expression = gamma * B

    # The final answer given by the LLM to be checked
    llm_final_answer = "<<<C>>>"

    # --- Verification Logic ---

    # 1. Extract the chosen option letter from the LLM's answer string
    try:
        chosen_option_letter = llm_final_answer.strip().split('<<<')[1].split('>>>')[0]
    except (IndexError, AttributeError):
        return "Invalid answer format. The answer must be in the format '<<<X>>>'."

    if chosen_option_letter not in options:
        return f"The chosen option '{chosen_option_letter}' is not one of the valid options (A, B, C, D)."

    # 2. Get the symbolic expression corresponding to the chosen option
    chosen_expression = options[chosen_option_letter]

    # 3. Compare the chosen expression with the correct physical expression
    #    sympy.simplify() is used to ensure the comparison is robust.
    if sympy.simplify(chosen_expression - correct_physical_expression) == 0:
        return "Correct"
    else:
        # Find which option letter corresponds to the correct physical answer
        correct_option_letter = None
        for letter, expression in options.items():
            if sympy.simplify(expression - correct_physical_expression) == 0:
                correct_option_letter = letter
                break
        
        reason = (f"The answer is incorrect. "
                  f"The correct oscillation frequency is the Larmor frequency, which is ω = γB. "
                  f"This corresponds to option {correct_option_letter} ({options[correct_option_letter]}). "
                  f"The provided answer chose option {chosen_option_letter}, which corresponds to the expression {chosen_expression}.")
        return reason

# Run the check and print the result
result = check_larmor_precession_answer()
print(result)