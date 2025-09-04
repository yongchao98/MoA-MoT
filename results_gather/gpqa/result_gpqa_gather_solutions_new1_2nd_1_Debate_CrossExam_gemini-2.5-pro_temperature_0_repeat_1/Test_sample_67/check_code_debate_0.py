import sympy

def check_larmor_frequency_answer():
    """
    Checks the correctness of the LLM's answer for the Larmor frequency problem.

    The problem asks for the oscillation frequency of the z-component of a magnetic
    moment that was initially aligned with the Z-axis after a magnetic field is
    suddenly switched to the Y-axis. The correct frequency is the Larmor frequency, ω = γB.
    """
    # 1. Define symbolic variables for a rigorous check
    gamma, B = sympy.symbols('gamma B', real=True, positive=True)

    # 2. State the correct physical formula for the oscillation frequency
    # This is the Larmor frequency.
    correct_frequency_formula = gamma * B

    # 3. Define the given multiple-choice options
    options = {
        'A': sympy.sqrt(2) * gamma * B,
        'B': gamma * B / 2,
        'C': gamma * B / sympy.sqrt(2),
        'D': gamma * B
    }

    # 4. The final answer provided by the LLM to be checked
    llm_answer_str = "<<<D>>>"

    # --- Verification Logic ---

    # Find the letter corresponding to the correct physical formula
    correct_option_letter = None
    for letter, formula in options.items():
        # Use sympy.simplify to check for mathematical equivalence
        if sympy.simplify(formula - correct_frequency_formula) == 0:
            correct_option_letter = letter
            break
    
    if correct_option_letter is None:
        # This case should not be reached if the correct answer is among the options
        return "Error: The correct physical answer (gamma*B) was not found in the options list."

    # Extract the letter from the LLM's answer string
    try:
        llm_choice_letter = llm_answer_str.strip().split('<<<')[1].split('>>>')[0]
    except IndexError:
        return f"The provided answer string '{llm_answer_str}' is not in the expected '<<<X>>>' format."

    # Compare the LLM's choice with the physically correct option
    if llm_choice_letter == correct_option_letter:
        return "Correct"
    else:
        reason = (
            f"The answer is incorrect.\n"
            f"Reason: The physical principle governing this system is Larmor precession. The angular frequency of oscillation, ω, is the Larmor frequency, given by the formula ω = γB.\n"
            f"The correct option is '{correct_option_letter}', as its expression '{options[correct_option_letter]}' matches the physical formula.\n"
            f"The provided answer was '{llm_choice_letter}', which corresponds to the expression '{options[llm_choice_letter]}'. This is physically incorrect."
        )
        return reason

# Execute the check and print the result
result = check_larmor_frequency_answer()
print(result)