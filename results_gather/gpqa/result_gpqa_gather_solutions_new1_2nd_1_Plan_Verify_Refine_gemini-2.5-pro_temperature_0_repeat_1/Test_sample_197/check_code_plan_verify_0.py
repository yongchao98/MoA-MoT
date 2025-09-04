import math

def check_correctness():
    """
    This function checks the correctness of the provided answer to the chemistry problem.
    It recalculates the percentage of the dithiocyanato cobalt(II) complex and compares it
    to the value of the chosen option.
    """
    
    # --- Problem Data ---
    # Stability constants
    beta1 = 9
    beta2 = 40
    beta3 = 63
    beta4 = 16
    
    # Ligand concentration
    L = 0.1  # [SCN-] in M

    # --- Options as presented in the final answer's context ---
    # A) 42.3%, B) 25.6%, C) 16.9%, D) 38.1%
    options = {
        'A': 42.3,
        'B': 25.6,
        'C': 16.9,
        'D': 38.1
    }
    
    # The final answer provided by the LLM
    llm_answer_choice = 'C'

    # --- Step 1: Recalculate the theoretical percentage ---
    # The formula for the mole fraction (alpha_n) of a complex [ML_n] is:
    # alpha_n = (beta_n * [L]^n) / (1 + sum(beta_i * [L]^i for i=1 to max_n))
    
    # Numerator for the [Co(SCN)2] species (n=2)
    numerator = beta2 * (L**2)
    
    # Denominator is the sum of terms for all species: Co^2+ (term=1), [Co(SCN)]+, [Co(SCN)2], etc.
    denominator = 1 + (beta1 * L**1) + (beta2 * L**2) + (beta3 * L**3) + (beta4 * L**4)
    
    # Check for division by zero, although unlikely in this context
    if denominator == 0:
        return "Calculation error: The denominator is zero."
        
    # Calculate the mole fraction and convert to percentage
    alpha2 = numerator / denominator
    calculated_percentage = alpha2 * 100

    # --- Step 2: Check the correctness of the LLM's answer ---
    
    # Constraint 1: The chosen option must be valid.
    if llm_answer_choice not in options:
        return f"Incorrect. The final answer '{llm_answer_choice}' is not a valid option (A, B, C, or D)."

    # Get the percentage value corresponding to the LLM's chosen option
    chosen_option_value = options[llm_answer_choice]

    # Constraint 2: The value of the chosen option must match the calculated value.
    # We use math.isclose() to handle potential floating-point inaccuracies.
    # A tolerance of 0.05% is reasonable for this type of problem.
    if math.isclose(calculated_percentage, chosen_option_value, abs_tol=0.05):
        return "Correct"
    else:
        # Find which option *should* have been chosen
        correct_option_choice = None
        for option, value in options.items():
            if math.isclose(calculated_percentage, value, abs_tol=0.05):
                correct_option_choice = option
                break
        
        if correct_option_choice:
            return (f"Incorrect. The calculation yields a percentage of approximately {calculated_percentage:.2f}%. "
                    f"This corresponds to option {correct_option_choice} ({options[correct_option_choice]}%). "
                    f"The provided answer selected option {llm_answer_choice} ({chosen_option_value}%).")
        else:
            return (f"Incorrect. The calculated percentage is approximately {calculated_percentage:.2f}%. "
                    f"This value does not match any of the provided options. "
                    f"The provided answer selected option {llm_answer_choice} ({chosen_option_value}%).")

# Run the check
result = check_correctness()
print(result)