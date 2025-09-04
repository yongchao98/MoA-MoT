import math

def check_answer_correctness():
    """
    Checks the correctness of the provided answer for the radiation problem.

    The problem asks for two quantities for an oscillating spheroidal charge distribution:
    1. The fraction of maximum power radiated at an angle theta = 30 degrees.
    2. A possible form of the power function f(lambda, theta), specifically its dependence on lambda.

    The provided answer is C, which corresponds to:
    - Fraction: 1/4
    - Wavelength dependence: proportional to lambda^(-4)
    """

    # --- Define the expected values from the answer (Option C) ---
    expected_fraction = 1/4
    expected_lambda_power = -4

    # --- Part 1: Verify the fraction of power ---
    # The standard model for this system is an oscillating electric dipole.
    # The power per unit solid angle, dP/dΩ, is proportional to sin^2(theta).
    # P(theta) ∝ sin^2(theta)

    # The maximum power 'A' occurs when sin^2(theta) is maximum.
    # The maximum value of sin(theta) is 1, which occurs at theta = 90 degrees.
    # So, A is proportional to sin^2(90°).
    theta_max_power_rad = math.radians(90)
    angular_factor_max = math.sin(theta_max_power_rad)**2  # This is 1.0

    # The power at the given angle theta = 30 degrees is proportional to sin^2(30°).
    theta_given_rad = math.radians(30)
    angular_factor_given = math.sin(theta_given_rad)**2  # This is (0.5)^2 = 0.25

    # The fraction of the maximum power is the ratio of the angular factors.
    # The proportionality constants cancel out.
    calculated_fraction = angular_factor_given / angular_factor_max

    # Check if the calculated fraction matches the expected fraction from the answer.
    if not math.isclose(calculated_fraction, expected_fraction):
        return (f"Incorrect: The fraction of power is wrong. "
                f"The power is proportional to sin^2(theta). The maximum power is at theta=90 deg (sin^2(90)=1). "
                f"At theta=30 deg, the power is proportional to sin^2(30) = (1/2)^2 = 0.25. "
                f"Therefore, the fraction should be 0.25 / 1 = 0.25. The answer states the fraction is {expected_fraction}, "
                f"but the calculation gives {calculated_fraction}.")

    # --- Part 2: Verify the wavelength dependence ---
    # For electric dipole radiation, the total radiated power is proportional to omega^4.
    # Since omega = 2*pi*c / lambda, the power is proportional to (1/lambda)^4 = lambda^-4.
    # This is the dominant radiation term, so it's the most "possible form" for f.
    model_lambda_power = -4

    if expected_lambda_power != model_lambda_power:
        return (f"Incorrect: The wavelength dependence is wrong. "
                f"For electric dipole radiation, power is proportional to lambda^{model_lambda_power}. "
                f"The answer suggests a dependence of lambda^{expected_lambda_power}, which does not match the standard physical model.")

    # If both checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result
result = check_answer_correctness()
print(result)