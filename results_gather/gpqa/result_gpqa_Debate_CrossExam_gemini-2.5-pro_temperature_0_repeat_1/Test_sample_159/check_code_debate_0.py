import numpy as np
from scipy.special import jn_zeros

def check_diffraction_answer():
    """
    Checks the correctness of the answer for the angular distance between the first two minima
    in Fraunhofer diffraction from a circular aperture.
    """
    # The problem simplifies to finding the angular distance between the first two minima
    # of the Airy diffraction pattern from a circular aperture of radius 'a'.

    # The angular position of the n-th minimum is given by:
    # theta_n ≈ z_n * lambda / (2 * pi * a)
    # where z_n are the zeros of the first-order Bessel function J1(x).

    # The angular distance is Delta_theta = theta_2 - theta_1.
    # Delta_theta = (z_2 - z_1) * lambda / (2 * pi * a)
    # We need to calculate the coefficient C = (z_2 - z_1) / (2 * pi).

    try:
        # Get the first two non-trivial zeros of J1(x) using SciPy for high precision.
        # The first argument is the order of the Bessel function (1).
        # The second argument is the number of zeros to find (2).
        zeros_j1 = jn_zeros(1, 2)
        z1 = zeros_j1[0]
        z2 = zeros_j1[1]
    except (ImportError, ModuleNotFoundError):
        # Fallback if scipy is not installed, using high-precision hardcoded values.
        print("SciPy not found. Using hardcoded values for Bessel function zeros.")
        z1 = 3.8317059702075127
        z2 = 7.015586669815618

    # Calculate the theoretical coefficient.
    theoretical_coefficient = (z2 - z1) / (2 * np.pi)

    # The options provided in the question correspond to these coefficients:
    options = {
        "A": 1.220,
        "B": 0.610,
        "C": 0.500,
        "D": 0.506
    }

    # The provided answer is D.
    provided_answer_key = "D"
    provided_answer_value = options[provided_answer_key]

    # The logic in the provided answer is sound:
    # 1. N-sided polygon with constant apothem 'a' -> circle of radius 'a' as N -> infinity.
    # 2. Diffraction from a circular aperture produces an Airy pattern.
    # 3. Minima are determined by the zeros of the J1 Bessel function.
    # 4. The formula for angular distance is correctly derived.
    # The final step is to check if the calculation leads to the chosen option.

    # We find which option's coefficient is closest to the theoretical value.
    # This accounts for potential rounding in the question's options.
    best_option_key = min(options, key=lambda k: abs(options[k] - theoretical_coefficient))

    if best_option_key == provided_answer_key:
        return "Correct"
    else:
        # This block would execute if the wrong option was chosen.
        return (
            f"The answer is incorrect because it selects the wrong option based on the calculation.\n"
            f"The theoretical coefficient is calculated as (z2 - z1) / (2 * pi) ≈ {theoretical_coefficient:.5f}.\n"
            f"The provided options are A={options['A']}, B={options['B']}, C={options['C']}, D={options['D']}.\n"
            f"The distances from the theoretical value are:\n"
            f"|A - calc| = {abs(options['A'] - theoretical_coefficient):.5f}\n"
            f"|B - calc| = {abs(options['B'] - theoretical_coefficient):.5f}\n"
            f"|C - calc| = {abs(options['C'] - theoretical_coefficient):.5f}\n"
            f"|D - calc| = {abs(options['D'] - theoretical_coefficient):.5f}\n"
            f"The closest option is '{best_option_key}' with a value of {options[best_option_key]}, but the provided answer was '{provided_answer_key}'."
        )

# Run the check
result = check_diffraction_answer()
print(result)