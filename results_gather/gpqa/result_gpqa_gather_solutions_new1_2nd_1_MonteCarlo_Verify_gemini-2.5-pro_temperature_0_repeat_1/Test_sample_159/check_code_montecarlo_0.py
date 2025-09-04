import numpy as np
from scipy.special import jn_zeros

def check_diffraction_answer():
    """
    Checks the correctness of the answer for the diffraction problem.

    The problem simplifies to finding the angular distance between the first two
    minima of a Fraunhofer diffraction pattern from a circular aperture of radius 'a'.
    The angular distance is given by Δθ = ((z₂ - z₁) / (2π)) * (λ/a), where z₁ and z₂
    are the first two non-trivial zeros of the Bessel function J₁(x).

    The code calculates the coefficient C = (z₂ - z₁) / (2π) and compares it
    with the coefficient from the provided answer.
    """
    # The final answer given by the LLM is 'D', which corresponds to the coefficient 0.506
    # from the options list provided in the final analysis:
    # A) 1.220, B) 0.610, C) 0.500, D) 0.506
    llm_answer_coefficient = 0.506
    llm_answer_option = 'D'

    # Get the first two non-trivial zeros of the first-order Bessel function J₁(x)
    try:
        z1, z2 = jn_zeros(1, 2)
    except ImportError:
        return "Error: The 'scipy' library is required. Please install it using 'pip install scipy'."
    except Exception as e:
        return f"An unexpected error occurred: {e}"

    # Calculate the theoretical coefficient for the angular distance
    calculated_coefficient = (z2 - z1) / (2 * np.pi)

    # Check if the calculated coefficient matches the one from the answer
    # Use a tolerance for floating-point comparison
    tolerance = 0.001
    if abs(calculated_coefficient - llm_answer_coefficient) < tolerance:
        return "Correct"
    else:
        # Provide a detailed reason for the incorrectness
        reason = (
            f"The provided answer '{llm_answer_option}' corresponds to a coefficient of {llm_answer_coefficient}, "
            f"but the calculated correct coefficient is approximately {calculated_coefficient:.4f}.\n"
            f"The calculation is based on Δθ = ((z₂ - z₁) / (2π)) * (λ/a).\n"
            f"The first two zeros of J₁(x) are z₁≈{z1:.4f} and z₂≈{z2:.4f}.\n"
            f"Therefore, the coefficient is ({z2:.4f} - {z1:.4f}) / (2π) ≈ {calculated_coefficient:.4f}.\n"
        )
        
        # Check for common mistakes (distractors)
        coeff_first_min = z1 / (2 * np.pi)
        if abs(coeff_first_min - 0.610) < tolerance:
            reason += (
                f"The coefficient for the position of the *first* minimum is ≈{coeff_first_min:.4f}, "
                f"which corresponds to option B (0.610). The question asks for the distance between the first two minima."
            )
        return reason

# Execute the check and print the result
print(check_diffraction_answer())