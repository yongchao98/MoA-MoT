import re
import math

def check_correctness_of_answer():
    """
    This function checks the correctness of the provided answer for the chemistry problem.
    It recalculates the percentage of the dithiocyanato cobalt(II) complex and compares it
    to the value corresponding to the selected option.
    """
    
    # --- Problem Data ---
    # Stability constants
    beta1 = 9
    beta2 = 40
    beta3 = 63
    beta4 = 16
    
    # Ligand concentration
    scn_conc = 0.1  # [SCN-] in M

    # Multiple choice options as defined in the final analysis block
    options = {
        "A": 38.1,
        "B": 16.9,
        "C": 42.3,
        "D": 25.6
    }

    # The final answer provided by the LLM to be checked
    llm_final_answer = "<<<B>>>"

    # --- Calculation ---
    # The formula for the mole fraction (alpha_n) of a complex [ML_n] is:
    # alpha_n = (beta_n * [L]^n) / (1 + sum(beta_i * [L]^i for i=1 to max_n))
    # We are calculating the percentage for n=2.

    # Numerator for the [Co(SCN)2] species
    numerator = beta2 * (scn_conc ** 2)

    # Denominator is the sum of terms for all species (Co^2+, [Co(SCN)]+, [Co(SCN)2], etc.)
    # The term for the free ion Co^2+ is 1.
    denominator = 1 + \
                  (beta1 * (scn_conc ** 1)) + \
                  (beta2 * (scn_conc ** 2)) + \
                  (beta3 * (scn_conc ** 3)) + \
                  (beta4 * (scn_conc ** 4))

    # Calculate the mole fraction (alpha_2)
    if denominator == 0:
        return "Error: Denominator in calculation is zero."
        
    alpha_2 = numerator / denominator

    # Convert to percentage
    calculated_percentage = alpha_2 * 100

    # --- Verification ---
    # Extract the letter from the LLM's final answer
    match = re.search(r'<<<([A-D])>>>', llm_final_answer)
    if not match:
        return f"Invalid answer format. Expected '<<<A>>>', '<<<B>>>', '<<<C>>>', or '<<<D>>>', but got {llm_final_answer}"

    answer_letter = match.group(1)
    
    # Get the numerical value corresponding to the chosen option
    expected_percentage = options.get(answer_letter)
    if expected_percentage is None:
        return f"Invalid option letter '{answer_letter}' found in the answer."

    # Compare the calculated result with the expected answer, allowing for a small tolerance
    # for rounding differences. A tolerance of 0.1% is reasonable given the options.
    if math.isclose(calculated_percentage, expected_percentage, abs_tol=0.1):
        return "Correct"
    else:
        return (f"Incorrect. The calculation shows the percentage of the dithiocyanato cobalt(II) complex "
                f"should be approximately {calculated_percentage:.2f}%. The provided answer is '{answer_letter}', "
                f"which corresponds to {expected_percentage}%. The calculated value does not match the answer.")

# Run the check and print the result
print(check_correctness_of_answer())