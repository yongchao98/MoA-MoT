import math

def check_answer_correctness():
    """
    This function checks the correctness of the provided answer to the physics question.
    It models the oscillating spheroid as an electric dipole and calculates the
    expected fraction of power and the wavelength dependence.
    """
    
    # The final answer provided by the LLM is 'C'.
    # Let's define what option 'C' represents.
    # Option C: fraction = 1/4, lambda_dependence = lambda^-4
    expected_fraction = 1/4
    expected_lambda_power = -4

    # --- Part 1: Calculate the fraction of maximum power at theta = 30 degrees ---
    
    # For an electric dipole, the radiated power per unit solid angle is proportional to sin^2(theta).
    # P(theta) ∝ sin^2(theta)
    
    # The maximum power (A) occurs at theta = 90 degrees, where sin^2(90) = 1.
    # A ∝ sin^2(90) = 1
    
    # The power at theta = 30 degrees is proportional to sin^2(30).
    angle_in_degrees = 30
    angle_in_radians = math.radians(angle_in_degrees)
    
    # The fraction of the maximum power is the ratio P(30)/P(90) = sin^2(30)/sin^2(90)
    calculated_fraction = math.sin(angle_in_radians)**2
    
    # --- Part 2: Determine the wavelength dependence ---
    
    # The power radiated by an oscillating electric dipole is proportional to the
    # fourth power of the angular frequency (ω^4).
    # Since ω = 2*pi*c/λ, power is proportional to (1/λ)^4, which is λ^-4.
    calculated_lambda_power = -4
    
    # --- Part 3: Compare calculated values with the provided answer (Option C) ---
    
    # Check if the calculated fraction matches the expected fraction from option C.
    # We use math.isclose for robust floating-point comparison.
    fraction_matches = math.isclose(calculated_fraction, expected_fraction)
    
    # Check if the calculated lambda power matches the expected power from option C.
    lambda_power_matches = (calculated_lambda_power == expected_lambda_power)
    
    if fraction_matches and lambda_power_matches:
        return "Correct"
    else:
        reasons = []
        if not fraction_matches:
            reasons.append(
                f"The fraction of maximum power is incorrect. "
                f"Based on the electric dipole model (Power ∝ sin^2(θ)), the fraction at 30 degrees is sin^2(30°) = (1/2)^2 = 0.25. "
                f"The answer's fraction is {expected_fraction}, but the calculation yields {calculated_fraction:.4f}."
            )
        if not lambda_power_matches:
            reasons.append(
                f"The wavelength dependence is incorrect. "
                f"Based on the electric dipole model (Power ∝ ω^4), the dependence on wavelength should be λ^{calculated_lambda_power}. "
                f"The answer's dependence is λ^{expected_lambda_power}."
            )
        return "Incorrect. " + " ".join(reasons)

# Run the check and print the result
result = check_answer_correctness()
print(result)