import scipy.special as sp
import numpy as np
import math

def check_diffraction_answer():
    """
    Checks the correctness of the answer for the angular distance between
    the first two diffraction minima of a circular aperture.
    """
    # The problem requires calculating the coefficient C in the expression:
    # Delta_theta = C * (lambda / a)
    # The theoretical value of C is (z2 - z1) / (2 * pi), where z1 and z2
    # are the first two non-trivial zeros of the Bessel function J1(x).

    try:
        # Get the first two non-trivial zeros of J1(x) using SciPy for high precision.
        zeros = sp.jn_zeros(1, 2)
        z1 = zeros[0]
        z2 = zeros[1]
    except ImportError:
        return "SciPy library not found. Cannot perform the check. Please install it using 'pip install scipy'."
    except Exception as e:
        return f"An error occurred while getting Bessel zeros: {e}"

    # Calculate the theoretical coefficient.
    theoretical_coefficient = (z2 - z1) / (2 * np.pi)

    # The LLM's final answer is 'C', which corresponds to the value 0.506.
    # The options as listed in the final analysis are:
    # A) 0.500, B) 1.220, C) 0.506, D) 0.610
    llm_answer_value = 0.506

    # The options are given to 3 decimal places. We check if the theoretical value
    # is closest to the selected option. A tolerance of 0.001 is appropriate.
    if math.isclose(theoretical_coefficient, llm_answer_value, abs_tol=0.001):
        return "Correct"
    else:
        # Provide a detailed reason if the check fails.
        # Calculate coefficients for other incorrect but plausible options.
        first_min_pos_coeff = z1 / (2 * np.pi)  # Corresponds to option D
        
        reason = (
            f"The answer is incorrect.\n"
            f"The LLM's chosen answer corresponds to a coefficient of {llm_answer_value}.\n"
            f"The correct theoretical coefficient is (z2 - z1) / (2 * pi) ≈ {theoretical_coefficient:.5f}.\n"
            f"The difference between the theoretical value and the answer's value is {abs(theoretical_coefficient - llm_answer_value):.5f}, which is larger than the acceptable tolerance for the given options.\n\n"
            f"For context:\n"
            f"- The coefficient for the position of the *first* minimum is z1 / (2 * pi) ≈ {first_min_pos_coeff:.5f}, which corresponds to option D (0.610).\n"
            f"- The correct answer, {theoretical_coefficient:.5f}, is numerically closest to 0.506 (Option C), so the LLM's choice was correct, but the check failed, indicating a potential issue with the tolerance or the provided answer value."
        )
        return reason

# Run the check
result = check_diffraction_answer()
print(result)