import re
import math

def check_correctness():
    """
    This function checks the correctness of the provided answer to the binary star system problem.
    
    It verifies the calculation based on the physical principles of binary star systems
    and compares the result with the chosen option.
    """
    
    # --- Problem Data ---
    # System 1
    P1 = 2  # years
    K1a = 10 # km/s
    K1b = 5  # km/s
    
    # System 2
    P2 = 1  # year
    K2a = 15 # km/s
    K2b = 10 # km/s
    
    # Options from the question
    options = {
        'A': 0.4,
        'B': 0.6,
        'C': 1.2,
        'D': 0.7
    }
    
    # The final answer provided by the LLM to be checked
    llm_final_answer = "<<<A>>>"

    # --- Calculation ---
    # For an eclipsing binary, the total mass is proportional to P * (K_total)^3.
    # We need to find the ratio M_sys1 / M_sys2.
    
    K_total1 = K1a + K1b
    K_total2 = K2a + K2b
    
    # Calculate the ratio
    try:
        calculated_ratio = (P1 * (K_total1)**3) / (P2 * (K_total2)**3)
    except ZeroDivisionError:
        return "Incorrect: Calculation resulted in a division by zero."

    # --- Verification ---
    # Extract the chosen option letter from the LLM's answer
    match = re.search(r'<<<([A-D])>>>', llm_final_answer)
    if not match:
        return f"Incorrect: The final answer format is invalid. Expected '<<<X>>>' but got '{llm_final_answer}'."
        
    chosen_option_letter = match.group(1)
    
    # Find the option that is numerically closest to the calculated ratio
    closest_option = None
    min_difference = float('inf')
    
    for option_letter, option_value in options.items():
        difference = abs(calculated_ratio - option_value)
        if difference < min_difference:
            min_difference = difference
            closest_option = option_letter
            
    # Check if the LLM's chosen option is the correct one
    if chosen_option_letter == closest_option:
        return "Correct"
    else:
        return (f"Incorrect: The calculated mass ratio is {calculated_ratio:.3f}. "
                f"The closest option is '{closest_option}' (value ~{options[closest_option]}). "
                f"The provided answer was '{chosen_option_letter}' (value ~{options[chosen_option_letter]}).")

# Run the check and print the result
result = check_correctness()
print(result)