import numpy as np
from scipy.special import jn_zeros
import re

def check_diffraction_answer():
    """
    Checks the correctness of the provided answer for the diffraction problem.

    The problem asks for the angular distance between the first two minima
    of a Fraunhofer diffraction pattern from a circular aperture of radius 'a'.

    The formula for the angular position of the m-th minimum is:
    θ_m = (z_m * λ) / (2 * π * a)
    where z_m is the m-th non-trivial zero of the Bessel function J₁(x).

    The angular distance between the first two minima is:
    Δθ = θ₂ - θ₁ = (z₂ - z₁) * λ / (2 * π * a)

    The code calculates the coefficient C = (z₂ - z₁) / (2 * π) and compares it
    with the coefficient from the chosen answer.
    """
    # The final answer provided by the LLM
    llm_answer_text = "<<<B>>>"

    # Define the options from the question
    options = {
        "A": 0.610,
        "B": 0.506,
        "C": 1.220,
        "D": 0.500
    }

    # --- Step 1: Perform the correct physical calculation ---

    # Get the first two non-trivial zeros of the Bessel function of the first kind, order 1 (J1).
    # jn_zeros(1, 2) returns an array of the first 2 zeros of J1.
    try:
        z1, z2 = jn_zeros(1, 2)
    except ImportError:
        return "Could not verify the answer because 'scipy' is not installed. Please install it using 'pip install scipy'."
    
    # Calculate the coefficient for the angular distance between the first two minima.
    # C = (z₂ - z₁) / (2 * π)
    calculated_coefficient = (z2 - z1) / (2 * np.pi)

    # --- Step 2: Parse the LLM's chosen answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return f"Could not parse the final answer from the text: '{llm_answer_text}'"
    
    chosen_option_key = match.group(1)
    
    if chosen_option_key not in options:
        return f"The chosen option '{chosen_option_key}' is not a valid option (A, B, C, or D)."

    chosen_coefficient = options[chosen_option_key]

    # --- Step 3: Verify the correctness of the choice ---

    # Check if the chosen coefficient is the best match for the calculated one.
    # We use np.isclose for robust floating-point comparison.
    # The tolerance is set based on the precision of the options.
    if not np.isclose(calculated_coefficient, chosen_coefficient, atol=0.001):
        # Find the best matching option
        best_match_key = min(options, key=lambda k: abs(options[k] - calculated_coefficient))
        best_match_value = options[best_match_key]
        
        return (f"Incorrect. The provided answer chose option {chosen_option_key} with a coefficient of {chosen_coefficient}. "
                f"However, the correct calculation yields a coefficient of approximately {calculated_coefficient:.5f}. "
                f"The best match among the options is option {best_match_key} with a coefficient of {best_match_value}.")

    # --- Step 4: Verify the analysis of distractors (optional but good practice) ---
    
    # Check distractor A: position of the first minimum
    theta_1_coeff = z1 / (2 * np.pi)
    if not np.isclose(theta_1_coeff, options["A"], atol=0.001):
        return (f"The provided answer is correct, but its analysis of distractors might be incomplete. "
                f"The coefficient for the first minimum is {theta_1_coeff:.5f}, which should match option A (0.610).")

    # Check distractor C: angular diameter of the central maximum
    diameter_coeff = 2 * z1 / (2 * np.pi)
    if not np.isclose(diameter_coeff, options["C"], atol=0.001):
        return (f"The provided answer is correct, but its analysis of distractors might be incomplete. "
                f"The coefficient for the central maximum's diameter is {diameter_coeff:.5f}, which should match option C (1.220).")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_diffraction_answer()
print(result)