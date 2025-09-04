import math

def check_physics_answer():
    """
    Checks the correctness of the answer based on the principles of electric dipole radiation.
    """
    
    # --- Part 1: Verify the fraction based on angular dependence ---
    # The radiated power per unit solid angle is proportional to sin^2(theta).
    # The maximum power A is at theta = 90 degrees, where sin^2(90) = 1.
    # We need to find the fraction of A at theta = 30 degrees.
    
    theta_deg = 30
    theta_rad = math.radians(theta_deg)
    
    # The fraction is sin^2(30) / sin^2(90) = sin^2(30) / 1
    # math.sin(theta_rad)**2 gives sin^2(30)
    calculated_fraction = math.sin(theta_rad)**2
    
    # The expected fraction from the physical principle is 1/4.
    expected_fraction = 1/4
    
    # Use math.isclose for robust floating-point comparison.
    if not math.isclose(calculated_fraction, expected_fraction):
        return (f"Constraint check failed: The fraction of maximum power at 30 degrees is incorrect. "
                f"Based on sin^2(30), the fraction should be {expected_fraction}, but the calculation yields {calculated_fraction}.")

    # --- Part 2: Verify the wavelength dependence ---
    # The radiated power is proportional to omega^4 (angular frequency).
    # Angular frequency omega is inversely proportional to wavelength lambda (omega ~ 1/lambda).
    # Therefore, Power is proportional to (1/lambda)^4 = lambda^(-4).
    
    # The expected exponent for lambda is -4.
    expected_exponent = -4
    
    # --- Part 3: Check the provided answer against the derived values ---
    
    # The final answer provided by the LLM is 'A'.
    # Let's define the options from the question.
    # Options are in the format: (fraction, lambda_exponent)
    options = {
        "A": (1/4, -4),
        "B": (3/4, -6),
        "C": (1/2, -4),
        "D": (1/4, -3)
    }
    
    llm_answer_key = "A"
    
    # Get the values from the chosen option
    chosen_fraction, chosen_exponent = options[llm_answer_key]
    
    # Check if the fraction in the chosen option matches the derived fraction
    if not math.isclose(chosen_fraction, expected_fraction):
        return (f"Incorrect Answer: The fraction in the chosen answer '{llm_answer_key}' is {chosen_fraction}, "
                f"which does not satisfy the physical constraint. The correct fraction is {expected_fraction}.")
                
    # Check if the exponent in the chosen option matches the derived exponent
    if chosen_exponent != expected_exponent:
        return (f"Incorrect Answer: The wavelength dependence exponent in the chosen answer '{llm_answer_key}' is {chosen_exponent}, "
                f"which does not satisfy the physical constraint. The correct exponent is {expected_exponent}.")
                
    # If both parts of the chosen answer satisfy the constraints, the answer is correct.
    return "Correct"

# Execute the check and print the result
result = check_physics_answer()
print(result)