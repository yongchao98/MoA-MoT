import math
import re

def check_answer_correctness():
    """
    Checks the correctness of the LLM's answer based on physics principles.
    
    The problem involves electric dipole radiation, where:
    1. Radiated power per unit solid angle is proportional to sin^2(theta).
    2. Radiated power is proportional to omega^4, which means it's proportional to lambda^-4.
    """
    
    # Define the problem's options
    # Format: {Option_Key: (Fraction, Lambda_Power)}
    options = {
        'A': (1/2, -4),
        'B': (1/4, -3),
        'C': (3/4, -6),
        'D': (1/4, -4)
    }
    
    # The final answer provided by the LLM
    llm_final_answer_key = 'D'
    
    # Get the values from the chosen option
    chosen_fraction, chosen_lambda_power = options[llm_final_answer_key]
    
    # --- Verification Step 1: Calculate the correct power fraction ---
    # The fraction of maximum power at theta = 30 degrees is sin^2(30) / sin^2(90).
    # sin(90 deg) = 1, so the fraction is just sin^2(30 deg).
    theta_deg = 30
    theta_rad = math.radians(theta_deg)
    correct_fraction = math.sin(theta_rad)**2
    
    # --- Verification Step 2: Determine the correct wavelength dependence ---
    # For electric dipole radiation, power is proportional to omega^4.
    # Since omega = 2*pi*c/lambda, power is proportional to (1/lambda)^4 = lambda^-4.
    correct_lambda_power = -4
    
    # --- Compare and Conclude ---
    # Check if the fraction from the chosen option is correct
    fraction_is_correct = math.isclose(chosen_fraction, correct_fraction, rel_tol=1e-9)
    
    # Check if the lambda power from the chosen option is correct
    lambda_power_is_correct = (chosen_lambda_power == correct_lambda_power)
    
    if fraction_is_correct and lambda_power_is_correct:
        return "Correct"
    else:
        error_messages = []
        if not fraction_is_correct:
            error_messages.append(
                f"The power fraction is incorrect. The chosen option '{llm_final_answer_key}' has a fraction of {chosen_fraction}, "
                f"but the correct physical value is sin^2(30°) = {correct_fraction:.4f} (1/4)."
            )
        if not lambda_power_is_correct:
            error_messages.append(
                f"The wavelength dependence is incorrect. The chosen option '{llm_final_answer_key}' has a dependence of λ^({chosen_lambda_power}), "
                f"but the correct dependence for electric dipole radiation is λ^({correct_lambda_power})."
            )
        return "Incorrect: " + " ".join(error_messages)

# Run the check
result = check_answer_correctness()
print(result)