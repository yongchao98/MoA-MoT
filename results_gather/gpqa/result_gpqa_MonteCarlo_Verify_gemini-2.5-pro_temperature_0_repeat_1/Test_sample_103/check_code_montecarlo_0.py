import math

def check_exoplanet_period_ratio():
    """
    Checks the correctness of the LLM's answer for the exoplanet period ratio problem.
    """
    # --- Given data from the problem ---
    # Wavelength shift for planet #1's star
    delta_lambda_1 = 5.0  # in miliangstroms
    # Wavelength shift for planet #2's star
    delta_lambda_2 = 7.0  # in miliangstroms

    # --- Theoretical Derivation ---
    # The question asks for the ratio of the orbital periods: P2 / P1.
    # The key physical relationships are:
    # 1. The star's velocity amplitude K is proportional to the wavelength shift Δλ.
    #    Therefore, K1 / K2 = Δλ1 / Δλ2.
    # 2. For a given star mass and planet mass, the velocity amplitude K is
    #    proportional to the period P to the power of -1/3 (K ∝ P^(-1/3)).
    #    Therefore, K1 / K2 = (P1/P2)^(-1/3) = (P2/P1)^(1/3).
    #
    # Combining these gives: Δλ1 / Δλ2 = (P2 / P1)^(1/3).
    # To solve for the period ratio P2 / P1, we cube both sides:
    # P2 / P1 = (Δλ1 / Δλ2)^3

    # --- Calculation ---
    try:
        calculated_ratio = (delta_lambda_1 / delta_lambda_2) ** 3
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Verification ---
    # The LLM's answer is 'C'.
    llm_answer_option = 'C'
    options = {'A': 0.85, 'B': 1.40, 'C': 0.36, 'D': 1.96}

    if llm_answer_option not in options:
        return f"Incorrect. The provided answer '{llm_answer_option}' is not a valid option."

    # Find which of the given options is numerically closest to our calculated result.
    closest_option = min(options, key=lambda k: abs(options[k] - calculated_ratio))

    # Check if the LLM's chosen option is the one closest to the correct calculation.
    if llm_answer_option == closest_option:
        # The answer is correct. The calculated value is ~0.3644, which is closest to 0.36.
        return "Correct"
    else:
        # The answer is incorrect.
        return (f"Incorrect. The derivation P2/P1 = (Δλ1/Δλ2)^3 is correct. "
                f"Plugging in the values gives (5/7)^3 ≈ {calculated_ratio:.4f}. "
                f"This value is closest to option {closest_option} ({options[closest_option]}), "
                f"but the provided answer was option {llm_answer_option} ({options[llm_answer_option]}).")

# Run the check
result = check_exoplanet_period_ratio()
print(result)