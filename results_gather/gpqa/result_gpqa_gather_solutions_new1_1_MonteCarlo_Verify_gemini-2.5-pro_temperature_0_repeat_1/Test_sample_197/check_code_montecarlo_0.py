import math

def check_chemistry_calculation():
    """
    This function verifies the calculation for the percentage of the 
    dithiocyanato cobalt(II) complex based on the provided problem data.
    """
    
    # --- Problem Constraints and Data ---
    # Total cobalt concentration c(Co) = 10^-2 M (This is not needed for the percentage calculation)
    # Thiocyanate concentration [SCN-]
    scn_conc = 0.1  # M
    
    # Stability constants
    beta1 = 9
    beta2 = 40
    beta3 = 63
    beta4 = 16
    
    # Multiple choice options
    options = {
        'A': 38.1,
        'B': 25.6,
        'C': 42.3,
        'D': 16.9
    }
    
    # The final answer provided by the LLM to be checked
    llm_final_answer_letter = 'D'

    # --- Calculation based on chemical principles ---
    # The formula for the mole fraction (alpha_n) of a complex [ML_n] is:
    # alpha_n = (beta_n * [L]^n) / (1 + sum(beta_i * [L]^i for i=1 to N))
    # We need to find the percentage for n=2, which is alpha_2 * 100.

    try:
        # Calculate the numerator for the n=2 species
        numerator = beta2 * (scn_conc ** 2)

        # Calculate the denominator, which is the sum of terms for all species (n=0 to 4)
        # The term for n=0 (free Co2+) is 1.
        denominator = 1 + (beta1 * scn_conc**1) + (beta2 * scn_conc**2) + (beta3 * scn_conc**3) + (beta4 * scn_conc**4)

        if denominator == 0:
            return "Calculation Error: Denominator is zero."

        # Calculate the mole fraction (alpha_2)
        alpha_2 = numerator / denominator

        # Convert to percentage
        calculated_percentage = alpha_2 * 100

    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Verification ---
    # Get the value corresponding to the LLM's answer
    llm_answer_value = options.get(llm_final_answer_letter)
    
    if llm_answer_value is None:
        return f"The provided answer '{llm_final_answer_letter}' is not a valid option."

    # Check if the calculated percentage matches the LLM's chosen option's value.
    # We use math.isclose for robust floating-point comparison. A relative tolerance of 1% is suitable here.
    if math.isclose(calculated_percentage, llm_answer_value, rel_tol=0.01):
        return "Correct"
    else:
        # Find which option the calculation actually supports
        correct_option = 'None'
        for letter, value in options.items():
            if math.isclose(calculated_percentage, value, rel_tol=0.01):
                correct_option = letter
                break
        
        return (f"Incorrect. The calculation yields a percentage of approximately {calculated_percentage:.2f}%. "
                f"This value corresponds to option {correct_option} ({options.get(correct_option, 'N/A')}%). "
                f"The provided answer was {llm_final_answer_letter} ({llm_answer_value}%).")

# Run the check and print the result
result = check_chemistry_calculation()
print(result)