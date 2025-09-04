import numpy as np
from scipy.special import jn_zeros

def check_diffraction_answer():
    """
    Checks the correctness of the answer for the angular distance between the first two minima
    of a circular aperture diffraction pattern.
    """
    # The LLM's answer is B, which corresponds to a coefficient of 0.506.
    llm_answer_coefficient = 0.506

    # The problem states that for an N-sided polygon with N -> infinity,
    # the aperture becomes a circle with radius 'a' (the apothem).
    # The angular position of the nth minimum in the far-field (Fraunhofer) diffraction
    # pattern for a circular aperture is given by:
    # sin(theta_n) = z_n * lambda / (2 * pi * a)
    # where z_n are the zeros of the first-order Bessel function of the first kind, J1(x).

    # Using the small angle approximation (sin(theta) ≈ theta):
    # theta_n = z_n * lambda / (2 * pi * a)

    # The angular distance between the first two minima is:
    # delta_theta = theta_2 - theta_1
    # delta_theta = (z_2 - z_1) / (2 * pi) * (lambda / a)

    # We need to calculate the coefficient C = (z_2 - z_1) / (2 * pi).

    try:
        # Get the first two non-zero zeros of the first-order Bessel function J1(x).
        # jn_zeros(1, 2) computes the first 2 zeros of J1.
        zeros = jn_zeros(1, 2)
        z1 = zeros[0]
        z2 = zeros[1]
    except ImportError:
        # Fallback if scipy is not installed, using high-precision values.
        # This makes the script runnable even without scipy.
        print("Warning: SciPy library not found. Using hardcoded values for Bessel function zeros.")
        z1 = 3.831705970207512
        z2 = 7.015586669815619

    # Calculate the theoretical coefficient.
    calculated_coefficient = (z2 - z1) / (2 * np.pi)

    # The options are given to three decimal places, so a tolerance of 0.001 is appropriate.
    tolerance = 0.001

    if abs(calculated_coefficient - llm_answer_coefficient) < tolerance:
        return "Correct"
    else:
        reason = (
            f"The answer is incorrect.\n"
            f"The theoretical coefficient for λ/a is (z2 - z1) / (2 * pi).\n"
            f"Using the first two zeros of J1(x), z1 ≈ {z1:.4f} and z2 ≈ {z2:.4f}.\n"
            f"The calculated coefficient is ({z2:.4f} - {z1:.4f}) / (2 * pi) ≈ {calculated_coefficient:.4f}.\n"
            f"The coefficient from the chosen answer is {llm_answer_coefficient}.\n"
            f"The difference |{calculated_coefficient:.4f} - {llm_answer_coefficient}| is {abs(calculated_coefficient - llm_answer_coefficient):.4f}, "
            f"which is not within the tolerance of {tolerance}."
        )
        return reason

# Run the check
result = check_diffraction_answer()
print(result)