import math

def check_correctness_of_answer():
    """
    This function checks the correctness of the LLM's answer to the chemistry problem.
    It calculates the percentage of the dithiocyanato cobalt(II) complex based on the provided data
    and compares it to the option selected by the LLM.
    """

    # --- Problem Data ---
    # Given values from the question
    # Concentration of the thiocyanate ligand [SCN-] in M
    L = 0.1
    # Overall stability constants
    beta1 = 9
    beta2 = 40
    beta3 = 63
    beta4 = 16

    # Multiple choice options from the question
    options = {
        'A': 25.6,
        'B': 38.1,
        'C': 42.3,
        'D': 16.9
    }

    # The final answer provided by the LLM to be checked
    llm_answer_key = 'D'

    # --- Calculation ---
    # The percentage of a specific complex [MLn] is given by the formula:
    # % = ( (beta_n * [L]^n) / (1 + sum(beta_i * [L]^i)) ) * 100
    # For the dithiocyanato cobalt(II) complex, [Co(SCN)2], n=2.

    # Calculate the numerator term for the [Co(SCN)2] complex
    numerator = beta2 * (L**2)

    # Calculate the denominator, which is the sum of terms for all species
    # (Co^2+, [Co(SCN)]+, [Co(SCN)2], [Co(SCN)3]-, [Co(SCN)4]2-)
    denominator = 1 + (beta1 * L**1) + (beta2 * L**2) + (beta3 * L**3) + (beta4 * L**4)

    # Avoid division by zero, though unlikely in this context
    if denominator == 0:
        return "Calculation Error: Denominator is zero."

    # Calculate the mole fraction (alpha_2) and then the percentage
    calculated_percentage = (numerator / denominator) * 100

    # --- Verification ---
    # Get the percentage value corresponding to the LLM's chosen option
    llm_answer_value = options.get(llm_answer_key)

    if llm_answer_value is None:
        return f"Invalid answer key '{llm_answer_key}' provided. Must be one of {list(options.keys())}."

    # Check if the calculated percentage is close to the LLM's answer value.
    # A tolerance is used to account for rounding in the options.
    # The calculated value is ~16.916%. The option is 16.9%. The difference is ~0.016.
    # An absolute tolerance of 0.05 is sufficient to confirm they match.
    if math.isclose(calculated_percentage, llm_answer_value, abs_tol=0.05):
        return "Correct"
    else:
        # Find the closest option to the calculated value for a more informative error message
        closest_option_key = min(options, key=lambda k: abs(options[k] - calculated_percentage))
        
        reason = (f"Incorrect. The provided answer is {llm_answer_key} ({llm_answer_value}%), "
                  f"but the calculation based on the problem's data yields {calculated_percentage:.3f}%. "
                  f"The correct option should be {closest_option_key} ({options[closest_option_key]}%).")
        return reason

# Execute the check and print the result
result = check_correctness_of_answer()
print(result)