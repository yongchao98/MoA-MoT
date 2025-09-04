import numpy as np
from scipy.special import jn_zeros

def check_diffraction_answer():
    """
    This function checks the correctness of the provided answer to the physics problem.

    The problem asks for the angular distance between the first two minima of a
    Fraunhofer diffraction pattern from a circular aperture of radius 'a'.
    An N-sided polygon with N->infinity and constant apothem 'a' becomes a circle of radius 'a'.

    The angular position of the n-th minimum is given by:
    θ_n = (z_n * λ) / (2 * π * a)
    where z_n is the n-th non-zero root of the Bessel function J_1(x).

    The angular distance is Δθ = θ_2 - θ_1.
    Δθ = [(z_2 - z_1) / (2 * π)] * (λ / a)

    The code calculates the coefficient C = (z_2 - z_1) / (2 * π) and compares it
    to the coefficient from the selected answer option 'A'.
    """
    try:
        # The final answer provided in the prompt is 'A', which corresponds to the coefficient 0.506.
        expected_coefficient = 0.506

        # Get the first two positive zeros of the Bessel function of the first kind, order 1 (J1).
        # jn_zeros(1, 2) returns an array containing the first two zeros.
        zeros = jn_zeros(1, 2)
        z1 = zeros[0]
        z2 = zeros[1]

        # Calculate the theoretical coefficient for the angular distance.
        calculated_coefficient = (z2 - z1) / (2 * np.pi)

        # Check if the calculated coefficient matches the one from the answer 'A'.
        # A tolerance of 0.001 is reasonable since the option has 3 decimal places.
        if np.isclose(calculated_coefficient, expected_coefficient, atol=0.001):
            return "Correct"
        else:
            reason = "Incorrect.\n"
            reason += f"The theoretical coefficient for the angular distance Δθ = C * (λ/a) is calculated as C = (z₂ - z₁) / (2π).\n"
            reason += f"Using the first two non-zero roots of the J₁ Bessel function, z₁ ≈ {z1:.5f} and z₂ ≈ {z2:.5f}.\n"
            reason += f"The calculated coefficient is C ≈ ({z2:.5f} - {z1:.5f}) / (2π) ≈ {calculated_coefficient:.5f}.\n"
            reason += f"The coefficient from the provided answer 'A' is {expected_coefficient}.\n"
            reason += f"The calculated value {calculated_coefficient:.5f} does not match the provided answer's coefficient {expected_coefficient} within an acceptable tolerance."
            return reason

    except ImportError:
        return "Error: The 'scipy' or 'numpy' library is required to run this check. Please install it."
    except Exception as e:
        return f"An unexpected error occurred: {e}"

# Execute the check and print the result.
result = check_diffraction_answer()
print(result)