import math

def check_interplanar_distance():
    """
    Checks the correctness of the answer for the interplanar distance in a rhombohedral crystal.
    """
    # --- Problem Constraints and Given Data ---
    # Crystal system: Rhombohedral
    a = 10.0  # Interatomic distance in Angstrom
    alpha_deg = 30.0  # Angle in degrees
    h, k, l = 1, 1, 1 # Miller indices for the (111) plane

    # Multiple choice options provided in the question
    options = {
        'A': 9.08,
        'B': 9.54,
        'C': 10.05,
        'D': 8.95
    }

    # The final answer to be checked
    provided_answer_letter = 'B'
    
    # --- Calculation ---
    # For the (111) plane in a rhombohedral system, the formula simplifies to:
    # 1/d^2 = 3 / [a^2 * (1 + 2*cos(alpha))]
    
    # Convert angle from degrees to radians for trigonometric functions
    alpha_rad = math.radians(alpha_deg)
    
    # Calculate cos(alpha)
    cos_alpha = math.cos(alpha_rad)
    
    # Calculate 1/d^2
    try:
        denominator = (a**2) * (1 + 2 * cos_alpha)
        if denominator == 0:
            return "Error: Division by zero in the formula. The parameters are likely invalid."
        
        inv_d_sq = 3 / denominator
        
        if inv_d_sq <= 0:
            return "Error: Cannot take the square root of a non-positive number. The parameters lead to an invalid physical distance."
            
        # Calculate d
        calculated_d = 1 / math.sqrt(inv_d_sq)
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Verification ---
    # Get the value corresponding to the provided answer
    if provided_answer_letter not in options:
        return f"Incorrect. The provided answer letter '{provided_answer_letter}' is not a valid option."
        
    answer_value = options[provided_answer_letter]
    
    # Set a tolerance for comparing floating-point numbers. 
    # 0.01 is reasonable given the precision of the options.
    tolerance = 0.01

    # Check if the calculated value matches the answer's value within the tolerance
    if abs(calculated_d - answer_value) < tolerance:
        return "Correct"
    else:
        # If incorrect, find the best matching option to provide a more detailed reason.
        min_diff = float('inf')
        best_match_letter = None
        for letter, value in options.items():
            diff = abs(calculated_d - value)
            if diff < min_diff:
                min_diff = diff
                best_match_letter = letter
        
        reason = (f"Incorrect. The provided answer is {provided_answer_letter} ({answer_value} Angstrom), "
                  f"but the calculated interplanar distance is approximately {calculated_d:.4f} Angstrom. "
                  f"This calculated value is closest to option {best_match_letter} ({options[best_match_letter]} Angstrom).")
        return reason

# Run the check
result = check_interplanar_distance()
print(result)