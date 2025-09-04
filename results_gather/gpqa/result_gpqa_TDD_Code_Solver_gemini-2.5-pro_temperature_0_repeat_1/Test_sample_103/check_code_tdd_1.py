import math

def check_exoplanet_period_ratio():
    """
    Checks the correctness of the LLM's answer for the exoplanet period ratio problem.
    """
    # --- Given Information ---
    # Periodic wavelength shift for star #1 (in miliangstroms)
    delta_lambda_1 = 5.0
    # Periodic wavelength shift for star #2 (in miliangstroms)
    delta_lambda_2 = 7.0

    # --- LLM's Answer ---
    llm_answer_text_value = 0.36
    llm_selected_option = 'B'
    options = {
        'A': 1.40,
        'B': 0.36,
        'C': 0.85,
        'D': 1.96
    }

    # --- Physics Derivation & Calculation ---
    # The star's radial velocity amplitude (K) is proportional to the wavelength shift (Δλ).
    # K ∝ Δλ
    #
    # For a circular orbit, the orbital period (T) is inversely proportional to the cube of the
    # radial velocity amplitude (K), assuming star and planet masses are constant between systems.
    # T ∝ 1 / K³
    #
    # Therefore, the ratio of the periods T2/T1 is:
    # T2 / T1 = (K1 / K2)³ = (Δλ₁ / Δλ₂)³
    try:
        calculated_ratio = (delta_lambda_1 / delta_lambda_2) ** 3
    except ZeroDivisionError:
        return "Error: Division by zero. delta_lambda_2 cannot be zero."
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Verification ---
    # 1. Check if the numerical value in the LLM's text is correct.
    # We use a relative tolerance to account for potential rounding in the answer.
    if not math.isclose(calculated_ratio, llm_answer_text_value, rel_tol=0.02):
        return (f"Incorrect numerical value. The calculated ratio T2/T1 is (5/7)^3 ≈ {calculated_ratio:.4f}. "
                f"The provided answer value is {llm_answer_text_value}, which is not within a 2% tolerance.")

    # 2. Check if the selected option ('B') is correct.
    if llm_selected_option not in options:
        return f"Invalid option '{llm_selected_option}' provided. Valid options are A, B, C, D."

    selected_option_value = options[llm_selected_option]
    if not math.isclose(calculated_ratio, selected_option_value, rel_tol=0.02):
        return (f"Incorrect option selected. The calculated ratio is ≈ {calculated_ratio:.4f}. "
                f"The selected option '{llm_selected_option}' corresponds to the value {selected_option_value}, "
                f"which is incorrect.")

    # 3. Check for consistency between the text value and the option value.
    if not math.isclose(llm_answer_text_value, selected_option_value):
        return (f"Inconsistent answer. The text states the ratio is {llm_answer_text_value}, "
                f"but the selected option '{llm_selected_option}' corresponds to {selected_option_value}.")

    return "Correct"

# Run the check
result = check_exoplanet_period_ratio()
print(result)