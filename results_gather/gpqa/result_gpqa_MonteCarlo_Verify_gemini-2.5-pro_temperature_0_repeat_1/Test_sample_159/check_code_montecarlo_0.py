import numpy as np
from scipy.special import jn_zeros

def check_diffraction_answer():
    """
    Checks the correctness of the answer for the diffraction problem.
    The core of the problem is to calculate the coefficient C = (x2 - x1) / (2 * pi),
    where x1 and x2 are the first two non-zero roots of the Bessel function J1(x).
    """
    # The provided answer is A, which corresponds to a coefficient of 0.506.
    # The provided reasoning and derivation are correct. This code verifies the numerical result.
    
    try:
        # Use scipy.special.jn_zeros to find the first 2 positive roots of J1(x).
        # jn_zeros(n, k) finds the first k positive roots of the Bessel function J_n.
        zeros = jn_zeros(1, 2)
        x1 = zeros[0]  # First non-zero root of J1(x)
        x2 = zeros[1]  # Second non-zero root of J1(x)
    except Exception as e:
        return f"An error occurred while finding Bessel function zeros using scipy: {e}"

    # Calculate the coefficient for the angular distance formula: Δθ = C * (λ / a)
    calculated_coefficient = (x2 - x1) / (2 * np.pi)

    # The options provided in the question.
    options = {'A': 0.506, 'B': 1.220, 'C': 0.500, 'D': 0.610}
    
    # The answer given by the LLM is 'A'.
    llm_answer_key = 'A'
    llm_answer_value = options[llm_answer_key]

    # Find the option that is numerically closest to our calculation.
    best_match_key = min(options, key=lambda k: abs(options[k] - calculated_coefficient))

    # Check if the LLM's answer is the best match.
    # A small tolerance is used to account for rounding in the question's options.
    if best_match_key == llm_answer_key:
        return "Correct"
    else:
        return (f"The answer is incorrect. The reasoning is correct, but the chosen option is not the most accurate.\n"
                f"Calculation Steps:\n"
                f"1. The first non-zero root of J1(x) is x1 ≈ {x1:.5f}.\n"
                f"2. The second non-zero root of J1(x) is x2 ≈ {x2:.5f}.\n"
                f"3. The coefficient is C = (x2 - x1) / (2 * pi) = ({x2:.5f} - {x1:.5f}) / (2 * {np.pi:.5f}) ≈ {calculated_coefficient:.5f}.\n"
                f"4. The calculated coefficient {calculated_coefficient:.5f} is closest to option {best_match_key} ({options[best_match_key]}), not option {llm_answer_key} ({llm_answer_value}).")

# Execute the check
result = check_diffraction_answer()
print(result)