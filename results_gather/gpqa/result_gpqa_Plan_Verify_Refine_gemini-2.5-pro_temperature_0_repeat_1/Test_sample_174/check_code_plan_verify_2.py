import numpy as np

def check_answer():
    """
    Checks the correctness of the provided answer based on the physics of electric dipole radiation.
    """
    # The provided answer is C
    chosen_option = 'C'
    options = {
        'A': {'fraction': 1/4, 'lambda_power': -3},
        'B': {'fraction': 1/2, 'lambda_power': -4},
        'C': {'fraction': 1/4, 'lambda_power': -4},
        'D': {'fraction': 3/4, 'lambda_power': -6}
    }

    # --- Part 1: Verify the fraction of power at 30 degrees ---

    # The radiated power per unit solid angle is proportional to sin^2(theta).
    # Let's define a function representing this angular dependence.
    def angular_power_factor(theta_degrees):
        theta_radians = np.deg2rad(theta_degrees)
        return np.sin(theta_radians)**2

    # The maximum power 'A' occurs when sin^2(theta) is maximum, which is at theta = 90 degrees.
    max_power_factor = angular_power_factor(90)

    # The power at the specified angle, theta = 30 degrees.
    power_at_30_deg_factor = angular_power_factor(30)

    # The fraction of the maximum power is the ratio of the power factors.
    calculated_fraction = power_at_30_deg_factor / max_power_factor

    # Get the fraction from the chosen option.
    expected_fraction = options[chosen_option]['fraction']

    # Check if the calculated fraction matches the one in the answer.
    # We use np.isclose for safe floating-point comparison.
    if not np.isclose(calculated_fraction, expected_fraction):
        return (f"Incorrect fraction: The calculated fraction of power at 30 degrees is "
                f"{calculated_fraction:.4f} (sin^2(30)/sin^2(90)), but the answer states it is "
                f"{expected_fraction}.")

    # --- Part 2: Verify the wavelength dependence ---

    # According to the Larmor formula for an oscillating dipole, the radiated power
    # is proportional to the fourth power of the angular frequency (omega^4).
    # Since omega = 2*pi*c / lambda, the power is proportional to (1/lambda)^4,
    # which is lambda^(-4).
    calculated_lambda_power = -4

    # Get the lambda power from the chosen option.
    expected_lambda_power = options[chosen_option]['lambda_power']

    # Check if the lambda power matches the one in the answer.
    if calculated_lambda_power != expected_lambda_power:
        return (f"Incorrect wavelength dependence: The power should be proportional to "
                f"lambda^({calculated_lambda_power}), but the answer states it is "
                f"lambda^({expected_lambda_power}).")

    # If both checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result.
print(check_answer())