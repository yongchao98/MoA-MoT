import math

def check_binary_star_mass_ratio():
    """
    Checks the correctness of the answer for the binary star mass ratio problem.
    """
    # --- Step 1: Define the data from the question ---
    # System 1
    P1 = 2.0  # years
    K1_sum = 10.0 + 5.0  # km/s

    # System 2
    P2 = 1.0  # year
    K2_sum = 15.0 + 10.0  # km/s

    # Options as defined in the question prompt
    options = {
        "A": 0.7,
        "B": 0.6,
        "C": 1.2,
        "D": 0.4
    }
    
    # The final answer given by the LLM being checked
    llm_final_answer_letter = "D"

    # --- Step 2: Perform the calculation based on physics principles ---
    # For eclipsing binaries, the mass ratio M1/M2 is proportional to (P1 * K1_sum^3) / (P2 * K2_sum^3)
    try:
        calculated_ratio = (P1 * (K1_sum**3)) / (P2 * (K2_sum**3))
    except ZeroDivisionError:
        return "Error: Division by zero in calculation. P2 or K2_sum cannot be zero."
    
    # --- Step 3: Verify the numerical result ---
    # The reasoning in the answer correctly simplifies the calculation to 54/125
    expected_ratio = 54 / 125
    if not math.isclose(calculated_ratio, expected_ratio, rel_tol=1e-9):
        return (f"Calculation Mismatch. The code calculated the ratio as {calculated_ratio}, "
                f"but the expected value from the fraction 54/125 is {expected_ratio}.")

    # --- Step 4: Find the closest option ---
    # Determine which option letter corresponds to the value closest to our calculated ratio.
    closest_option_letter = min(options, key=lambda k: abs(options[k] - calculated_ratio))
    
    # --- Step 5: Compare the LLM's answer with the verified correct answer ---
    if llm_final_answer_letter == closest_option_letter:
        return "Correct"
    else:
        return (f"Incorrect. The calculated mass ratio is {calculated_ratio:.3f}. "
                f"The closest option value is {options[closest_option_letter]}, which corresponds to option '{closest_option_letter}'. "
                f"The provided answer was '{llm_final_answer_letter}'.")

# Run the check
result = check_binary_star_mass_ratio()
print(result)