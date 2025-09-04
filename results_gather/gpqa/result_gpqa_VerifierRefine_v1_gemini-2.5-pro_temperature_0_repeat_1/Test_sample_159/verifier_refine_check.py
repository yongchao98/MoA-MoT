import numpy as np
from scipy.special import jn_zeros

def check_diffraction_answer():
    """
    Checks the correctness of the answer for the diffraction problem.

    The problem asks for the angular distance between the first two minima for
    Fraunhofer diffraction through a circular aperture of radius 'a'.

    The angular position of the m-th minimum is given by:
    theta_m = z_m * lambda / (2 * pi * a)
    where z_m is the m-th non-zero root of the Bessel function J1(x).

    The angular distance between the first two minima is:
    Delta_theta = theta_2 - theta_1 = (z_2 - z_1) * lambda / (2 * pi * a)

    This function calculates the coefficient C = (z_2 - z_1) / (2 * pi) and
    compares it with the value from option A.
    """
    try:
        # Get the first two positive (non-zero) roots of the Bessel function J1(x).
        # jn_zeros(n, nt) computes the first nt positive zeros of Jn(x).
        zeros = jn_zeros(1, 2)
        z1 = zeros[0]
        z2 = zeros[1]
    except (ImportError, ModuleNotFoundError):
        # Fallback if scipy is not installed, using known high-precision values.
        print("Warning: scipy.special not found. Using pre-computed values for Bessel zeros.")
        z1 = 3.8317059702075123
        z2 = 7.015586669815618

    # The coefficient for the angular distance between the first two minima
    calculated_coefficient = (z2 - z1) / (2 * np.pi)

    # The value from option A
    option_A_value = 0.506

    # Check if the calculated coefficient matches option A within a reasonable tolerance
    if np.isclose(calculated_coefficient, option_A_value, atol=1e-3):
        # The answer is correct. Let's double-check the reasoning for other options.
        # Option B: Position of the first minimum, coefficient = z1 / (2*pi)
        option_B_value = 0.610
        first_min_coeff = z1 / (2 * np.pi)
        if not np.isclose(first_min_coeff, option_B_value, atol=1e-3):
             return (f"The main answer is correct, but the implied reasoning for option B is flawed. "
                     f"The coefficient for the first minimum is {first_min_coeff:.4f}, not exactly 0.610.")

        # Option C: Angular width of the central maximum, coefficient = 2 * z1 / (2*pi) = z1/pi
        option_C_value = 1.220
        central_max_width_coeff = z1 / np.pi
        if not np.isclose(central_max_width_coeff, option_C_value, atol=1e-3):
            return (f"The main answer is correct, but the implied reasoning for option C is flawed. "
                    f"The coefficient for the central maximum width is {central_max_width_coeff:.4f}, not exactly 1.220.")

        return "Correct"
    else:
        return (f"Incorrect. The answer claims the coefficient is {option_A_value}, but the "
                f"correct calculation yields a different value.\n"
                f"The first two non-zero roots of J1(x) are z1 = {z1:.4f} and z2 = {z2:.4f}.\n"
                f"The coefficient for the angular distance is (z2 - z1) / (2 * pi).\n"
                f"Calculation: ({z2:.4f} - {z1:.4f}) / (2 * {np.pi:.4f}) = {z2 - z1:.4f} / {2 * np.pi:.4f} = {calculated_coefficient:.4f}.\n"
                f"This value ({calculated_coefficient:.4f}) does not match the provided answer's coefficient of {option_A_value}.")

# Run the check and print the result
result = check_diffraction_answer()
print(result)