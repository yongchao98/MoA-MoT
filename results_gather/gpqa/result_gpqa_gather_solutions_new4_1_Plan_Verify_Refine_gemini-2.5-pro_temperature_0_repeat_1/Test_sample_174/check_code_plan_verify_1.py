import math

def check_answer_correctness():
    """
    This function checks the correctness of the provided answer for the physics question.
    The question asks for two quantities:
    1. The fraction of maximum power radiated at an angle of 30 degrees.
    2. The dependence of the radiated power on the wavelength (lambda).

    The final answer provided is 'A', which corresponds to the pair (1/4, lambda^-4).
    This code verifies these two values based on the principles of electric dipole radiation.
    """

    # --- Part 1: Check the fraction of maximum power ---

    # The standard model for an oscillating charge distribution is an electric dipole.
    # The power radiated per unit solid angle (dP/dOmega) is proportional to sin^2(theta),
    # where theta is the angle from the oscillation axis.
    # dP/dOmega = K * sin^2(theta)

    # The maximum power 'A' occurs when sin^2(theta) is maximum, which is 1 (at theta=90 deg).
    # So, A is proportional to K.

    # We calculate the power at theta = 30 degrees.
    angle_deg = 30.0
    angle_rad = math.radians(angle_deg)

    # The fraction of maximum power is the ratio of power at 30 deg to the max power.
    # Fraction = (K * sin^2(30)) / (K * sin^2(90)) = sin^2(30) / 1
    calculated_fraction = math.sin(angle_rad)**2
    expected_fraction = 1.0 / 4.0

    # Check if the calculated fraction matches the answer's fraction.
    if not math.isclose(calculated_fraction, expected_fraction):
        return (f"Incorrect: The fraction of maximum power is wrong. "
                f"The standard electric dipole model predicts a power dependence of sin^2(theta). "
                f"The fraction at 30 degrees should be sin^2(30) = (1/2)^2 = 0.25. "
                f"The answer's fraction is {expected_fraction}, but the calculation gives {calculated_fraction:.4f}.")

    # --- Part 2: Check the wavelength dependence ---

    # For electric dipole radiation, the radiated power is proportional to the fourth
    # power of the angular frequency (omega^4).
    # Since omega is inversely proportional to wavelength (omega ~ 1/lambda),
    # the power is proportional to (1/lambda)^4, or lambda^-4.
    theoretical_lambda_exponent = -4
    expected_lambda_exponent = -4 # From option A: lambda^-4

    if theoretical_lambda_exponent != expected_lambda_exponent:
        return (f"Incorrect: The wavelength dependence is wrong. "
                f"The standard electric dipole model predicts a power dependence of lambda^{theoretical_lambda_exponent}. "
                f"The answer states a dependence of lambda^{expected_lambda_exponent}.")

    # If both parts are correct according to the standard model, the answer is correct.
    return "Correct"

# Execute the check and print the result.
result = check_answer_correctness()
print(result)