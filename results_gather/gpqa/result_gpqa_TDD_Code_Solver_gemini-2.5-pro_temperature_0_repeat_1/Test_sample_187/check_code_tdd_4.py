import math

def check_correctness():
    """
    This function checks the correctness of the given answer to the physics problem.
    
    Problem Details:
    - Crystal system: Rhombohedral
    - Interatomic distance (a): 10 Angstrom
    - Interaxial angles (alpha = beta = gamma): 30 degrees
    - Miller indices (h, k, l): (1, 1, 1)
    - Question: What is the interplanar distance d_111?

    LLM's Answer:
    - B) 9.54 Angstrom
    """

    # --- Define problem parameters and the available options ---
    a = 10.0  # Angstrom
    alpha_deg = 30.0
    h, k, l = 1, 1, 1
    
    options = {
        "A": 8.95,
        "B": 9.54,
        "C": 10.05,
        "D": 9.08
    }
    
    # The LLM's selected answer key.
    llm_answer_key = "B"

    # --- Calculation ---
    # The formula for interplanar spacing d_hkl in a rhombohedral system is:
    # 1/d^2 = [ (h^2+k^2+l^2)sin^2(α) + 2(hk+kl+lh)(cos^2(α)-cos(α)) ] /
    #         [ a^2 * (1 - 3cos^2(α) + 2cos^3(α)) ]
    #
    # A more numerically stable version, which will be used here, is:
    # 1/d^2 = [ (h²+k²+l²)(1+cosα) - 2(hk+kl+lh)cosα ] / [ a²(1-cosα)(1+2cosα) ]

    try:
        alpha_rad = math.radians(alpha_deg)
        cos_a = math.cos(alpha_rad)

        # Check for potential division by zero in the denominator's factors.
        # This corresponds to invalid crystal structures (alpha=0 or 120 deg).
        if math.isclose(1 - cos_a, 0) or math.isclose(1 + 2 * cos_a, 0):
            return "Calculation failed: Denominator in the formula is zero or near-zero."

        # Calculate the terms for the numerator
        sum_hkl_sq = h**2 + k**2 + l**2
        sum_hkl_prod = h*k + k*l + l*h

        # Calculate numerator and denominator of the 1/d^2 formula
        numerator = sum_hkl_sq * (1 + cos_a) - 2 * sum_hkl_prod * cos_a
        denominator = a**2 * (1 - cos_a) * (1 + 2 * cos_a)

        inv_d_squared = numerator / denominator
        
        # The result of 1/d^2 must be positive to take a real square root.
        if inv_d_squared <= 0:
            return f"Calculation failed: 1/d^2 is non-positive ({inv_d_squared:.4f}), so a real distance cannot be calculated."

        # Calculate the final interplanar distance
        calculated_d = 1 / math.sqrt(inv_d_squared)

    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Verification ---
    # To check a multiple-choice answer, we find which option is closest to our calculated value.
    # Then we check if this matches the LLM's chosen option.
    
    # Find the absolute difference between the calculated value and each option's value.
    errors = {key: abs(value - calculated_d) for key, value in options.items()}
    
    # Find the option key with the smallest error.
    best_match_key = min(errors, key=errors.get)

    # Compare the best matching option with the LLM's answer.
    if best_match_key == llm_answer_key:
        return "Correct"
    else:
        return (f"Incorrect. The calculated interplanar distance is {calculated_d:.4f} Angstrom. "
                f"This value is closest to option {best_match_key} ({options[best_match_key]} Angstrom), "
                f"but the provided answer was option {llm_answer_key} ({options[llm_answer_key]} Angstrom).")

# Run the check and print the result.
result = check_correctness()
print(result)