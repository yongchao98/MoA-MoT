import math

def check_answer():
    """
    Checks the correctness of the LLM's answer for the cobalt thiocyanate complex problem.
    """
    # Given values from the question
    SCN_conc = 0.1  # M
    beta1 = 9
    beta2 = 40
    beta3 = 63
    beta4 = 16

    # The formula for the mole fraction of the dithiocyanato complex (alpha_2) is:
    # alpha_2 = (beta_2 * [SCN-]^2) / (1 + beta_1*[SCN-]^1 + beta_2*[SCN-]^2 + beta_3*[SCN-]^3 + beta_4*[SCN-]^4)

    # Calculate the numerator (term for [Co(SCN)2])
    numerator = beta2 * (SCN_conc ** 2)

    # Calculate the denominator (sum of terms for all species)
    term0 = 1  # for Co(II)
    term1 = beta1 * (SCN_conc ** 1)
    term2 = beta2 * (SCN_conc ** 2)
    term3 = beta3 * (SCN_conc ** 3)
    term4 = beta4 * (SCN_conc ** 4)
    denominator = term0 + term1 + term2 + term3 + term4
    
    # Verify the intermediate calculations from the LLM's response
    expected_numerator = 0.4
    if not math.isclose(numerator, expected_numerator, rel_tol=1e-9):
        return f"Incorrect calculation: The numerator was calculated as {numerator}, but it should be {expected_numerator}."

    expected_denominator = 2.3646
    if not math.isclose(denominator, expected_denominator, rel_tol=1e-9):
        return f"Incorrect calculation: The denominator was calculated as {denominator}, but it should be {expected_denominator}."

    # Calculate the mole fraction and percentage
    if denominator == 0:
        return "Error: Denominator is zero, cannot perform division."
    
    alpha2 = numerator / denominator
    calculated_percentage = alpha2 * 100

    # The options provided in the final LLM response
    options = {
        'A': 38.1,
        'B': 25.6,
        'C': 42.3,
        'D': 16.9
    }

    # The final answer provided by the LLM
    llm_answer_letter = 'D'
    
    # Check if the letter is a valid option
    if llm_answer_letter not in options:
        return f"Invalid answer format: The final answer '{llm_answer_letter}' is not one of the options {list(options.keys())}."

    # Get the percentage value corresponding to the LLM's answer
    llm_answer_value = options[llm_answer_letter]

    # Check if the calculated percentage matches the value of the chosen option
    # Using a tolerance to account for floating point inaccuracies and rounding
    if math.isclose(calculated_percentage, llm_answer_value, rel_tol=1e-2, abs_tol=0.05):
        return "Correct"
    else:
        return (f"Incorrect. The calculated percentage is approximately {calculated_percentage:.1f}%. "
                f"The LLM's answer is '{llm_answer_letter}', which corresponds to {llm_answer_value}%. "
                f"The calculated value does not match the chosen option's value.")

# Run the check
result = check_answer()
print(result)