import math

def check_correctness_of_answer():
    """
    This function checks the correctness of the LLM's answer to the physics problem.
    It verifies the two components of the answer:
    1. The fraction of maximum power radiated at 30 degrees.
    2. The dependence of radiated power on the wavelength.
    """
    
    # --- Step 1: Define the problem parameters and the LLM's answer ---
    
    # The LLM's final answer choice is extracted from the provided text.
    llm_answer_choice = 'C'
    
    # Define the options provided in the question.
    options = {
        'A': {'fraction': 1/4, 'lambda_exp': -3},
        'B': {'fraction': 1/2, 'lambda_exp': -4},
        'C': {'fraction': 1/4, 'lambda_exp': -4},
        'D': {'fraction': 3/4, 'lambda_exp': -6}
    }
    
    # --- Step 2: Calculate the correct values based on physics principles ---
    
    # Part A: Calculate the fraction of maximum power.
    # The physical model is an oscillating electric dipole.
    # The power radiated per unit solid angle (dP/dOmega) is proportional to sin^2(theta).
    # The maximum power (A) occurs at theta = 90 degrees, where sin^2(90) = 1.
    # We need to find the power at theta = 30 degrees relative to this maximum.
    
    theta_deg = 30
    theta_rad = math.radians(theta_deg)
    
    # The power factor at 30 degrees is sin^2(30).
    power_factor_at_30 = math.sin(theta_rad)**2
    
    # The maximum power factor is sin^2(90) = 1.
    max_power_factor = 1.0
    
    # The fraction is the ratio of the power factors.
    calculated_fraction = power_factor_at_30 / max_power_factor
    
    # Part B: Determine the wavelength dependence.
    # For electric dipole radiation, the radiated power is proportional to the fourth
    # power of the angular frequency (omega^4).
    # Since omega is inversely proportional to wavelength lambda (omega ~ 1/lambda),
    # the power is proportional to (1/lambda)^4, which is lambda^(-4).
    calculated_lambda_exponent = -4
    
    # --- Step 3: Compare the calculated values with the LLM's chosen option ---
    
    if llm_answer_choice not in options:
        return f"Error: The provided answer choice '{llm_answer_choice}' is not a valid option (A, B, C, or D)."
        
    expected_values = options[llm_answer_choice]
    expected_fraction = expected_values['fraction']
    expected_lambda_exp = expected_values['lambda_exp']
    
    # Check if both parts of the answer are correct.
    is_fraction_correct = math.isclose(calculated_fraction, expected_fraction)
    is_lambda_exp_correct = (calculated_lambda_exponent == expected_lambda_exp)
    
    if is_fraction_correct and is_lambda_exp_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_fraction_correct:
            error_messages.append(
                f"The fraction of maximum power is incorrect. "
                f"The correct fraction, based on sin^2(30), is {calculated_fraction:.2f} (or 1/4). "
                f"The chosen option '{llm_answer_choice}' states the fraction is {expected_fraction}."
            )
        if not is_lambda_exp_correct:
            error_messages.append(
                f"The wavelength dependence is incorrect. "
                f"The correct dependence for electric dipole radiation is lambda^{calculated_lambda_exponent}. "
                f"The chosen option '{llm_answer_choice}' states the dependence is lambda^{expected_lambda_exp}."
            )
        return "\n".join(error_messages)

# Execute the check.
result = check_correctness_of_answer()
print(result)