import numpy as np
from scipy.special import jn_zeros

def check_diffraction_answer():
    """
    Checks the correctness of the answer for the diffraction problem.

    The problem asks for the angular distance between the first two minima for
    diffraction through a circular aperture of radius 'a'. The formula for this
    distance is Δθ = C * (λ/a), where the coefficient C = (z₂ - z₁) / (2π).
    z₁ and z₂ are the first two non-trivial zeros of the first-order Bessel
    function J₁(x).

    This code calculates the theoretical coefficient C and compares it to the
    coefficient from the selected answer 'B' (0.506).
    """

    # The options for the coefficient C are given in the question.
    options = {
        "A": 0.610,
        "B": 0.506,
        "C": 1.220,
        "D": 0.500
    }
    
    # The final answer provided by the LLM is 'B'.
    final_answer_key = "B"
    
    if final_answer_key not in options:
        return f"Invalid answer key '{final_answer_key}'. Valid keys are {list(options.keys())}."
        
    final_answer_value = options[final_answer_key]

    try:
        # Find the first two non-trivial zeros of J₁(x).
        # jn_zeros(n, nt) finds the first nt zeros of J_n(x).
        zeros = jn_zeros(1, 2)
        z1 = zeros[0]  # First zero, approx 3.8317
        z2 = zeros[1]  # Second zero, approx 7.0156

        # Calculate the theoretical coefficient for the angular distance Δθ = θ₂ - θ₁
        calculated_coefficient = (z2 - z1) / (2 * np.pi)

        # Find which option is numerically closest to the calculated value.
        # This is the correct way to select from multiple-choice options.
        closest_option_key = min(options, key=lambda k: abs(options[k] - calculated_coefficient))

        # Check if the given answer key matches the key of the closest option.
        if closest_option_key == final_answer_key:
            return "Correct"
        else:
            # Provide a detailed reason for the incorrectness.
            # Also check for common distractors.
            coeff_theta1 = z1 / (2 * np.pi) # Coefficient for the first minimum
            
            reason = (f"Incorrect. The provided answer is '{final_answer_key}' with a coefficient of {final_answer_value}. "
                      f"The theoretical coefficient is calculated as (z₂ - z₁) / (2π), which is approximately {calculated_coefficient:.5f}. "
                      f"The closest option to this value is '{closest_option_key}' with a coefficient of {options[closest_option_key]}.")
            
            if final_answer_value == round(coeff_theta1, 3):
                 reason += f"\nThe selected answer {final_answer_value} corresponds to the position of the first minimum (θ₁), not the distance between the first two minima (θ₂ - θ₁)."

            return reason

    except ImportError:
        return "Verification failed: The 'scipy' library is required but not installed. Please install it using 'pip install scipy'."
    except Exception as e:
        return f"An error occurred during verification: {e}"

# Execute the check and print the result.
result = check_diffraction_answer()
print(result)