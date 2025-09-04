import math

def check_rhombohedral_distance():
    """
    Checks the correctness of the calculated interplanar distance for a rhombohedral crystal.
    """
    # --- Given Parameters ---
    a = 10.0  # Lattice parameter in Angstrom
    alpha_deg = 30.0  # Interaxial angle in degrees
    h, k, l = 1, 1, 1 # Miller indices for the (111) plane

    # --- Options and Answer from the LLM response to be checked ---
    # The final response provides these options and selects 'A'.
    options = {
        "A": 9.54,
        "B": 10.05,
        "C": 9.08,
        "D": 8.95
    }
    llm_selected_option = "A"
    llm_answer_value = options[llm_selected_option]

    # --- Calculation ---
    # Using the simplified formula for the (111) plane, which is more robust:
    # d_111 = a * sqrt((1 + 2*cos(alpha)) / 3)
    try:
        # Convert angle from degrees to radians for Python's math functions
        alpha_rad = math.radians(alpha_deg)
        
        # Calculate the term inside the square root
        term_in_sqrt = (1 + 2 * math.cos(alpha_rad)) / 3
        
        if term_in_sqrt < 0:
            return "Calculation Error: The term inside the square root is negative."
            
        # Calculate the final interplanar distance
        calculated_d = a * math.sqrt(term_in_sqrt)

    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Verification ---
    # Check if the calculated value is close to the value from the selected option.
    # A tolerance is used for floating-point comparisons.
    tolerance = 0.01  # A tolerance of 0.01 Angstrom is reasonable.
    
    if abs(calculated_d - llm_answer_value) <= tolerance:
        # The numerical value is correct. Now, confirm it's the best match.
        best_match_letter = None
        min_difference = float('inf')
        for letter, value in options.items():
            difference = abs(calculated_d - value)
            if difference < min_difference:
                min_difference = difference
                best_match_letter = letter
        
        # Check if the LLM's chosen letter corresponds to the best-matching option
        if best_match_letter == llm_selected_option:
            return "Correct"
        else:
            return (f"The final answer letter '{llm_selected_option}' is incorrect. "
                    f"The calculated value is {calculated_d:.4f} Angstrom, which best matches "
                    f"option '{best_match_letter}' ({options[best_match_letter]} Angstrom), not option '{llm_selected_option}'.")
    else:
        return (f"The numerical value in the selected answer '{llm_selected_option}' ({llm_answer_value} Angstrom) is incorrect. "
                f"The calculated interplanar distance is {calculated_d:.4f} Angstrom.")

# Execute the check
result = check_rhombohedral_distance()
print(result)