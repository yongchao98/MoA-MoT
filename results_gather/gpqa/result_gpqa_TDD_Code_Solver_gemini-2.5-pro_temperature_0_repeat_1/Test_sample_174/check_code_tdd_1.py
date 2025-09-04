import math

def check_correctness():
    """
    This function checks the correctness of the provided answer for a physics problem
    concerning the radiation from an oscillating spheroid charge distribution.

    The question asks for two quantities:
    1. The fraction of maximum power radiated at an angle θ = 30°.
    2. The functional form of the power's dependence on wavelength λ.

    The provided answer (Option D) states:
    - Fraction = 1/4
    - Functional form ∝ λ⁻⁴
    """

    # --- Part 1: Verification of the power fraction ---

    # The physical model for an oscillating spheroid is an oscillating electric dipole.
    # The power radiated per unit solid angle (dP/dΩ) is proportional to sin²(θ).
    # The maximum power is radiated at θ = 90°, where sin²(90°) = 1.
    # The fraction of the maximum power at θ = 30° is therefore sin²(30°).
    
    try:
        theta_deg = 30.0
        theta_rad = math.radians(theta_deg)
        
        # Calculate the expected fraction based on physics principles.
        expected_fraction = math.sin(theta_rad)**2
        
        # Get the fraction from the provided answer (Option D).
        answer_fraction = 1/4
        
        # Check if the answer's fraction matches the expected value.
        if not math.isclose(expected_fraction, answer_fraction):
            return (f"Constraint not satisfied: The power fraction is incorrect. "
                    f"For an electric dipole, the fraction of maximum power at 30 degrees "
                    f"is sin²(30°) = (1/2)² = 0.25. The answer's fraction is {answer_fraction}, "
                    f"which does not match the calculated value of {expected_fraction:.2f}.")

    except Exception as e:
        return f"An error occurred during the fraction check: {e}"

    # --- Part 2: Verification of the wavelength dependence ---

    # According to Larmor's formula, the total power (P) radiated by an oscillating dipole
    # is proportional to the fourth power of the angular frequency (ω), i.e., P ∝ ω⁴.
    # Angular frequency (ω) is inversely proportional to wavelength (λ), i.e., ω ∝ λ⁻¹.
    # Therefore, the power P must be proportional to (λ⁻¹)⁴ = λ⁻⁴.
    
    try:
        # The exponent of λ in the answer's functional form (λ⁻⁴).
        answer_lambda_exponent = -4
        
        # The expected exponent based on physics principles.
        power_vs_omega_exponent = 4
        omega_vs_lambda_exponent = -1
        expected_lambda_exponent = power_vs_omega_exponent * omega_vs_lambda_exponent
        
        # Check if the answer's exponent matches the expected one.
        if answer_lambda_exponent != expected_lambda_exponent:
            return (f"Constraint not satisfied: The wavelength dependence is incorrect. "
                    f"Radiated power is proportional to ω⁴ and ω ∝ λ⁻¹, which means power should be "
                    f"proportional to λ⁻⁴. The answer suggests a dependence of λ^({answer_lambda_exponent}), "
                    f"which does not match the expected exponent of {expected_lambda_exponent}.")

    except Exception as e:
        return f"An error occurred during the wavelength dependence check: {e}"

    # If both parts of the answer are verified successfully, the answer is correct.
    return "Correct"

# Run the check.
result = check_correctness()
print(result)