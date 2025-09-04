import numpy as np

def check_interplanar_distance():
    """
    This function calculates the interplanar distance for a rhombohedral crystal
    and checks if the provided answer is correct.
    """
    # --- Problem Constraints & Given Values ---
    # Lattice parameter 'a' in Angstrom
    a = 10.0
    # Interaxial angle 'alpha' in degrees
    alpha_deg = 30.0
    # Miller indices (h, k, l) for the plane
    h, k, l = 1, 1, 1

    # --- The LLM's Answer ---
    llm_answer_option = 'B'
    llm_answer_value = 9.54

    # --- Calculation ---
    # Convert angle from degrees to radians for trigonometric functions
    alpha_rad = np.deg2rad(alpha_deg)
    cos_a = np.cos(alpha_rad)
    sin_a = np.sin(alpha_rad)

    # The formula for interplanar spacing in a rhombohedral system:
    # 1/d^2 = [ (h^2+k^2+l^2)sin^2(a) + 2(hk+kl+lh)(cos^2(a)-cos(a)) ] /
    #         [ a^2 * (1 - 3cos^2(a) + 2cos^3(a)) ]

    # Calculate the numerator of the formula
    term1_num = (h**2 + k**2 + l**2) * sin_a**2
    term2_num = 2 * (h*k + k*l + l*h) * (cos_a**2 - cos_a)
    numerator = term1_num + term2_num

    # Calculate the denominator of the formula
    # The volume term (1 - 3cos^2(a) + 2cos^3(a)) can also be written as
    # (1 - cos(a))^2 * (1 + 2cos(a)), which helps avoid floating point issues.
    denominator_factor = 1 - 3*cos_a**2 + 2*cos_a**3
    denominator = a**2 * denominator_factor

    # Ensure the denominator is valid (not zero or negative)
    if denominator <= 1e-9:
        print(f"Calculation error: The denominator is {denominator}, which is too close to zero. This can happen for invalid crystal parameters.")
        return

    # Calculate 1/d^2 and then d
    d_squared_inv = numerator / denominator
    calculated_d = 1.0 / np.sqrt(d_squared_inv)

    # --- Verification ---
    # Check if the calculated value is close to the LLM's answer value
    # A tolerance of 0.01 is reasonable given the precision of the options.
    tolerance = 0.01
    if abs(calculated_d - llm_answer_value) < tolerance:
        print("Correct")
    else:
        # If the answer is wrong, explain why.
        print(f"Incorrect. The provided answer is {llm_answer_value} Angstrom (Option {llm_answer_option}).")
        print(f"The correct interplanar distance, calculated using the standard formula, is {calculated_d:.4f} Angstrom.")
        
        # Find the best matching option from the list
        options = {'A': 10.05, 'B': 9.54, 'C': 8.95, 'D': 9.08}
        best_match_key = min(options.keys(), key=lambda key: abs(options[key] - calculated_d))
        best_match_value = options[best_match_key]
        
        print(f"The calculated value is closest to Option {best_match_key} ({best_match_value} Angstrom).")
        if best_match_key != llm_answer_option:
             print(f"The LLM chose the wrong option. It should have been {best_match_key}.")

# Run the check
check_interplanar_distance()