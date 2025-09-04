import numpy as np
try:
    # For a precise calculation, we use scipy's special function to find the zeros of the Bessel function.
    from scipy.special import jn_zeros
    scipy_available = True
except ImportError:
    # If scipy is not installed, we can use pre-calculated high-precision values.
    scipy_available = False

def check_correctness():
    """
    This function checks the correctness of the given answer by calculating the theoretical value.
    
    The problem simplifies to finding the angular distance between the first two minima
    of the Fraunhofer diffraction pattern from a circular aperture of radius 'a'.
    
    The angular position of the n-th minimum is given by:
    theta_n = z_n * lambda / (2 * pi * a)
    where z_n is the n-th zero of the J1 Bessel function.
    
    The angular distance between the first two minima is:
    delta_theta = theta_2 - theta_1 = (z_2 - z_1) * lambda / (2 * pi * a)
    
    The code calculates the coefficient C = (z_2 - z_1) / (2 * pi) and compares it
    to the coefficient from the given answer 'A', which is 0.506.
    """
    
    # The coefficient from the given answer 'A' is 0.506.
    given_answer_coefficient = 0.506

    if scipy_available:
        # Get the first two positive zeros of the Bessel function of the first kind, order 1 (J1).
        zeros = jn_zeros(1, 2)
        z1, z2 = zeros[0], zeros[1]
    else:
        # Fallback to known high-precision values if scipy is not available.
        z1 = 3.8317059702075123  # First zero of J1
        z2 = 7.015586669815612   # Second zero of J1

    # Calculate the theoretical coefficient for the term (lambda / a).
    calculated_coefficient = (z2 - z1) / (2 * np.pi)

    # Check if the calculated coefficient is close to the given answer's coefficient.
    # A tolerance of 1e-3 is suitable as the options are given to 3 decimal places.
    if np.isclose(calculated_coefficient, given_answer_coefficient, atol=1e-3):
        return "Correct"
    else:
        return (f"The answer is incorrect. The given answer corresponds to a coefficient of {given_answer_coefficient}, "
                f"but the theoretically calculated coefficient is {calculated_coefficient:.5f}.\n"
                f"Calculation: (z2 - z1) / (2*pi) = ({z2:.5f} - {z1:.5f}) / (2*pi) = {calculated_coefficient:.5f}.")

# The final output is the result of the check.
# print(check_correctness())