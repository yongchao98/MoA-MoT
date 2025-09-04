import numpy as np
import re

def check_answer_correctness():
    """
    Checks the correctness of the given answer (Option A) based on the physics
    of electric dipole radiation.
    """
    # The provided answer is Option A.
    # Let's define the parameters for this option.
    answer_fraction = 1/4
    answer_lambda_form = "λ^(-4)"

    # --- Part 1: Verify the fraction of power at 30 degrees ---

    # According to the theory of electric dipole radiation, the power radiated
    # per unit solid angle is proportional to sin^2(theta).
    # P(theta) ∝ sin^2(theta)
    # The maximum power 'A' occurs at theta = 90 degrees, where sin^2(90) = 1.
    # The power at theta = 30 degrees is proportional to sin^2(30).
    # The fraction of the maximum power is sin^2(30) / sin^2(90).

    theta_degrees = 30
    theta_radians = np.deg2rad(theta_degrees)
    
    # sin(30 degrees) = 1/2, so sin^2(30 degrees) = 1/4
    expected_fraction = np.sin(theta_radians)**2

    # Check if the answer's fraction matches the calculated physical value.
    # We use np.isclose for safe floating-point comparison.
    if not np.isclose(answer_fraction, expected_fraction):
        return (f"Incorrect: The fraction of power is wrong. "
                f"The physical model predicts a fraction of sin^2(30°) = {expected_fraction:.4f}. "
                f"The answer provides a fraction of {answer_fraction}.")

    # --- Part 2: Verify the wavelength dependence ---

    # The power radiated by an oscillating dipole is proportional to the fourth
    # power of the angular frequency (ω^4).
    # P ∝ ω^4
    # The relationship between angular frequency (ω) and wavelength (λ) is ω = 2πc/λ.
    # Therefore, the power is proportional to (1/λ)^4, which is λ^(-4).
    # The expected exponent for λ is -4.

    expected_exponent = -4
    
    # We parse the exponent from the string provided in the answer.
    try:
        # Use regular expression to find the integer in λ^(-d) or λ^(-d)
        match = re.search(r'λ\^\(?(-?\d+)\)?', answer_lambda_form)
        if match:
            answer_exponent = int(match.group(1))
        else:
            return f"Incorrect: Could not parse the exponent from the lambda dependence string '{answer_lambda_form}'."
    except (ValueError, IndexError):
        return f"Incorrect: Could not parse the exponent from the lambda dependence string '{answer_lambda_form}'."

    # Check if the answer's exponent matches the expected exponent.
    if answer_exponent != expected_exponent:
        return (f"Incorrect: The wavelength dependence is wrong. "
                f"The physical model predicts a power dependence of λ^({expected_exponent}). "
                f"The answer implies a dependence of λ^({answer_exponent}).")

    # If both checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
result = check_answer_correctness()
print(result)