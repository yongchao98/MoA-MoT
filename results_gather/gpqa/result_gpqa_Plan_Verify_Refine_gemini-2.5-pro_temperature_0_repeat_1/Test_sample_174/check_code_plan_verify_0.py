import numpy as np
import math

def check_dipole_radiation_answer():
    """
    Checks the correctness of the provided answer for the oscillating spheroid problem.
    The problem is modeled as electric dipole radiation.
    The answer to check is Option C: Fraction = 1/4, Wavelength dependence = lambda^(-4).
    """

    # --- Part 1: Verify the fraction of power at theta = 30 degrees ---

    # The power radiated per unit solid angle is proportional to sin^2(theta).
    # The maximum power (A) is at theta = 90 degrees, where sin^2(90) = 1.
    # The fraction of A at theta = 30 degrees is sin^2(30).
    
    theta_deg = 30.0
    # Convert degrees to radians for numpy's sin function
    theta_rad = np.deg2rad(theta_deg)
    
    # Calculate the theoretical fraction
    calculated_fraction = np.sin(theta_rad)**2
    
    # The fraction from the chosen answer (Option C) is 1/4
    answer_fraction = 1.0 / 4.0
    
    # Check if the calculated fraction matches the answer's fraction
    if not math.isclose(calculated_fraction, answer_fraction):
        return (f"Incorrect: The fraction of power at 30 degrees is wrong. "
                f"The power is proportional to sin^2(theta). At 30 degrees, the fraction of maximum power "
                f"should be sin^2(30) = (1/2)^2 = {calculated_fraction:.2f}. "
                f"The answer states the fraction is {answer_fraction}, which does not match.")

    # --- Part 2: Verify the wavelength dependence ---

    # The radiated power (P) is proportional to the fourth power of the angular frequency (omega^4).
    # P ~ omega^4
    # Angular frequency (omega) is inversely proportional to wavelength (lambda).
    # omega ~ 1/lambda
    # Therefore, the power P is proportional to (1/lambda)^4 = lambda^(-4).
    
    # The theoretical exponent for lambda is -4.
    theoretical_exponent = -4
    
    # The exponent from the chosen answer (Option C) is -4.
    answer_exponent = -4
    
    if theoretical_exponent != answer_exponent:
        return (f"Incorrect: The wavelength dependence is wrong. "
                f"Radiated power is proportional to omega^4, and since omega is proportional to 1/lambda, "
                f"the power should be proportional to lambda^({theoretical_exponent}). "
                f"The answer suggests a dependence of lambda^({answer_exponent}).")

    # --- Conclusion ---
    # If both checks pass, the answer is correct.
    return "Correct"

# Execute the check
result = check_dipole_radiation_answer()
print(result)