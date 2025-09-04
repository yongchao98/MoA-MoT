import math

def check_exoplanet_period_ratio():
    """
    This function checks the correctness of the provided answer for the exoplanet period ratio problem.
    """
    
    # --- Step 1: Define problem constants and the provided answer ---
    
    # Wavelength shift for planet #1's system in miliangstrom
    delta_lambda_1 = 5.0
    
    # Wavelength shift for planet #2's system in miliangstrom
    delta_lambda_2 = 7.0
    
    # The multiple-choice options provided in the question
    options = {
        "A": 1.40,
        "B": 0.36,
        "C": 0.85,
        "D": 1.96
    }
    
    # The answer selected by the LLM
    llm_selected_option = "B"

    # --- Step 2: Re-derive the physical relationship to verify the formula ---
    
    # The radial velocity amplitude of the star (K) is proportional to the observed wavelength shift (Δλ).
    # K ∝ Δλ
    
    # For a planet in a circular orbit, the star's velocity amplitude (K) is related to the system parameters by:
    # K ∝ m_p * M_star^(-2/3) * T^(-1/3)
    # where m_p is planet mass, M_star is star mass, and T is the orbital period.
    
    # The problem states that M_star and m_p are the same for both systems.
    # Therefore, these terms are constant and can be ignored in the proportionality.
    # This simplifies the relationship to:
    # K ∝ T^(-1/3)
    
    # Combining the two proportionalities:
    # Δλ ∝ K ∝ T^(-1/3)
    # This implies that Δλ * T^(1/3) is a constant for both systems.
    # So, Δλ1 * T1^(1/3) = Δλ2 * T2^(1/3)
    
    # We need to find the ratio T2 / T1. Let's rearrange the equation:
    # (T2 / T1)^(1/3) = Δλ1 / Δλ2
    # T2 / T1 = (Δλ1 / Δλ2)³
    
    # The formula T2/T1 = (Δλ1/Δλ2)³ used in the provided answer is correct.

    # --- Step 3: Calculate the expected value and compare with the selected answer ---
    
    try:
        # Calculate the ratio using the derived formula
        calculated_ratio = (delta_lambda_1 / delta_lambda_2) ** 3
    except ZeroDivisionError:
        return "Error: delta_lambda_2 cannot be zero."
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # Get the numerical value of the LLM's chosen answer
    llm_answer_value = options[llm_selected_option]
    
    # Check if the calculated result is close to the value of the chosen option.
    # The options are rounded, so we check for proximity, not exact equality.
    # A relative tolerance of 5% is reasonable for multiple-choice questions.
    if not math.isclose(calculated_ratio, llm_answer_value, rel_tol=0.05):
        # If it's not close, find the option that is actually the closest.
        closest_option = min(options, key=lambda k: abs(options[k] - calculated_ratio))
        return (f"Incorrect. The calculated ratio T2/T1 is {calculated_ratio:.4f}. "
                f"This value is closest to option {closest_option} ({options[closest_option]}), "
                f"but the provided answer was option {llm_selected_option} ({llm_answer_value}).")

    # --- Step 4: Final verification ---
    # The derivation is correct, the formula is correct, the calculation is correct,
    # and the result matches the selected option.
    return "Correct"

# Run the check and print the result
result = check_exoplanet_period_ratio()
print(result)