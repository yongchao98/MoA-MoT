import math

def check_answer_correctness():
    """
    This function checks the correctness of the answer to the physics question.
    It verifies the two components of the answer: the power fraction and the
    wavelength dependence, based on the principles of electric dipole radiation.
    """

    # --- Part 1: Verify the fraction of radiated power ---

    # The power radiated per unit solid angle for an electric dipole is proportional to sin²(θ).
    # The maximum power (A) is at θ=90°, where sin²(90°)=1.
    # The fraction of power at θ=30° is sin²(30°).
    angle_degrees = 30
    angle_radians = math.radians(angle_degrees)
    
    # Calculate the theoretical fraction
    theoretical_fraction = math.sin(angle_radians)**2
    
    # --- Part 2: Verify the wavelength dependence ---

    # The power radiated by an oscillating electric dipole is proportional to ω⁴.
    # Since ω is proportional to 1/λ, the power is proportional to (1/λ)⁴ = λ⁻⁴.
    # Therefore, the theoretical exponent of λ is -4.
    theoretical_lambda_exponent = -4

    # --- Part 3: Compare with the given options ---
    
    # The options provided in the question are:
    # A) 3/4, λ⁻⁶
    # B) 1/4, λ⁻³
    # C) 1/4, λ⁻⁴
    # D) 1/2, λ⁻⁴
    
    # The final answer from the LLM is 'C'.
    # Let's check if option 'C' matches our theoretical derivation.
    answer_c_fraction = 1/4
    answer_c_exponent = -4

    # Check if the fraction from option C is correct
    if not math.isclose(theoretical_fraction, answer_c_fraction):
        return (f"Incorrect: The fraction of power is wrong. "
                f"Based on the sin²(θ) law, the fraction at 30 degrees should be sin²(30°) = {theoretical_fraction:.2f}. "
                f"Option C states the fraction is {answer_c_fraction}, which is correct, but the check failed unexpectedly.")

    # Check if the wavelength dependence from option C is correct
    if theoretical_lambda_exponent != answer_c_exponent:
        return (f"Incorrect: The wavelength dependence is wrong. "
                f"Based on dipole radiation theory, power is proportional to λ⁻⁴ (exponent -4). "
                f"Option C states the exponent is {answer_c_exponent}, which is correct, but the check failed unexpectedly.")

    # If both parts of option C are consistent with the physical theory, the answer is correct.
    return "Correct"

# Run the check
result = check_answer_correctness()
print(result)