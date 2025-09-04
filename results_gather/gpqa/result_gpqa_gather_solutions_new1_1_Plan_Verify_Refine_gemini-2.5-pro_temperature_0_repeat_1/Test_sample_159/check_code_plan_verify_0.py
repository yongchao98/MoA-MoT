import numpy as np
from scipy.special import jn_zeros

def check_diffraction_answer():
    """
    Checks the correctness of the provided answer for the diffraction problem.
    The final answer from the LLM is 'A', which corresponds to the coefficient 0.506.
    """
    # --- Problem Parameters and LLM Answer ---
    # The options are A) 0.506, B) 0.610, C) 1.220, D) 0.500
    options = {'A': 0.506, 'B': 0.610, 'C': 1.220, 'D': 0.500}
    llm_answer_option = 'A'
    llm_answer_coefficient = options[llm_answer_option]

    # --- Theoretical Calculation ---
    # The core of the problem is to calculate the coefficient C in the expression:
    # Δθ = C * (λ / a)
    # where C = (z₂ - z₁) / (2 * π).
    # z₁ and z₂ are the first and second non-trivial zeros of the J1 Bessel function.

    try:
        # Get the first two non-trivial zeros of the Bessel function of the first kind, order 1 (J1).
        zeros = jn_zeros(1, 2)
        z1 = zeros[0]  # First zero
        z2 = zeros[1]  # Second zero
    except ImportError:
        return "Could not verify the answer because 'scipy' is not installed. Please install it using 'pip install scipy'."
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # Calculate the theoretical coefficient
    calculated_coefficient = (z2 - z1) / (2 * np.pi)

    # --- Verification ---
    # Find the option with the minimum absolute difference from the calculated value
    best_match_option = min(options, key=lambda k: abs(options[k] - calculated_coefficient))

    # Check if the LLM's chosen option is the best match
    if llm_answer_option == best_match_option:
        # The LLM correctly identified the best option among the choices.
        # Now, check if the reasoning is sound by comparing the calculated value.
        # A tolerance of 0.001 is reasonable given the precision of the options.
        if abs(calculated_coefficient - llm_answer_coefficient) < 0.001:
            return "Correct"
        else:
            # This case is unlikely if the option is the best match, but included for completeness.
            return (f"The chosen option {llm_answer_option} is the closest match, but the numerical value differs significantly. "
                    f"Calculated coefficient: {calculated_coefficient:.5f}, "
                    f"Option {llm_answer_option} value: {llm_answer_coefficient}.")
    else:
        return (f"Incorrect. The final answer is wrong. "
                f"The calculated coefficient for the angular distance is {calculated_coefficient:.5f} * (λ/a). "
                f"This value is closest to option {best_match_option} ({options[best_match_option]}). "
                f"The provided answer selected option {llm_answer_option} ({llm_answer_coefficient}), which is not the best match.")

# Execute the check and print the result
result = check_diffraction_answer()
print(result)