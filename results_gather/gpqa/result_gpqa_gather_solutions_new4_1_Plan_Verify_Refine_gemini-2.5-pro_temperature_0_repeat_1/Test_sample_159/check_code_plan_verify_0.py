import numpy as np
from scipy.special import jn_zeros

def check_diffraction_answer():
    """
    This function verifies the answer to the diffraction problem.

    The problem asks for the angular distance between the first two minima for a
    circular aperture of radius 'a'. The formula for this distance is:
    Δθ = (z₂ - z₁) / (2π) * (λ / a)
    where z₁ and z₂ are the first two non-trivial zeros of the J₁(x) Bessel function.

    The code calculates the coefficient C = (z₂ - z₁) / (2π) and checks if it
    matches the coefficient from the provided answer.
    """
    try:
        # The provided answer is 'B', which corresponds to a coefficient of 0.506.
        provided_answer_label = 'B'
        options = {
            "A": 0.610,
            "B": 0.506,
            "C": 0.500,
            "D": 1.220
        }
        
        # Get the first two non-trivial zeros of the J1 Bessel function.
        # jn_zeros(1, 2) returns the first 2 zeros of J1(x).
        zeros = jn_zeros(1, 2)
        z1 = zeros[0]  # approx 3.831706
        z2 = zeros[1]  # approx 7.015587

        # Calculate the correct coefficient for the angular distance.
        correct_coefficient = (z2 - z1) / (2 * np.pi)

        # Find the option that is numerically closest to the calculated correct value.
        closest_option_label = min(options.keys(), key=lambda k: abs(options[k] - correct_coefficient))

        # Check if the provided answer matches the closest correct option.
        if provided_answer_label == closest_option_label:
            return "Correct"
        else:
            # Explain why the answer is wrong.
            reason = (
                f"Incorrect. The provided answer is '{provided_answer_label}', which corresponds to a coefficient of {options[provided_answer_label]}.\n"
                f"The correct calculation for the coefficient is (z₂ - z₁) / (2π), where z₁ and z₂ are the first two zeros of the J₁ Bessel function.\n"
                f"z₁ ≈ {z1:.4f}, z₂ ≈ {z2:.4f}\n"
                f"Calculated coefficient = ({z2:.4f} - {z1:.4f}) / (2π) ≈ {correct_coefficient:.4f}.\n"
                f"This value is closest to option '{closest_option_label}' ({options[closest_option_label]}), not '{provided_answer_label}'."
            )
            return reason

    except ImportError:
        return "Could not perform check: The 'scipy' or 'numpy' library is not installed."
    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check
result = check_diffraction_answer()
print(result)