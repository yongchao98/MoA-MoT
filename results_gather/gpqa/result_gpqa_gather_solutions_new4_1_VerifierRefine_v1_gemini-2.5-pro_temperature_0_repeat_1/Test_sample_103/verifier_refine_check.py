import math

def check_answer():
    """
    Checks the correctness of the provided answer for the exoplanet orbital period problem.
    """
    # --- Problem Constraints & Given Data ---
    # Wavelength shift for planet #1's host star (in miliangstrom)
    delta_lambda_1 = 5
    # Wavelength shift for planet #2's host star (in miliangstrom)
    delta_lambda_2 = 7

    # The options as presented in the question
    options = {
        'A': 1.96,
        'B': 0.85,
        'C': 0.36,
        'D': 1.40
    }

    # The final answer provided by the LLM to be checked
    llm_answer_choice = 'C'
    
    # --- Derivation from Physical Principles ---
    # The radial velocity semi-amplitude (K) is proportional to the wavelength shift (Δλ).
    # K ∝ Δλ
    # The relationship between K and the orbital period (P) for the given constraints
    # (same star mass, same planet mass, assuming same inclination) is:
    # K ∝ P^(-1/3)
    # From this, we can derive the relationship between the periods and the wavelength shifts:
    # P ∝ K^(-3) ∝ (Δλ)^(-3)
    # Therefore, the ratio of the periods P₂/P₁ is:
    # P₂ / P₁ = (Δλ₂)⁻³ / (Δλ₁)⁻³ = (Δλ₁ / Δλ₂)³

    try:
        # Calculate the theoretical ratio P₂ / P₁
        calculated_ratio = (delta_lambda_1 / delta_lambda_2) ** 3
    except Exception as e:
        return f"An error occurred during the calculation: {e}"

    # --- Verification ---
    # Check if the LLM's chosen option exists
    if llm_answer_choice not in options:
        return f"The provided answer choice '{llm_answer_choice}' is not one of the valid options {list(options.keys())}."

    # Get the numerical value corresponding to the LLM's choice
    llm_answer_value = options[llm_answer_choice]

    # Compare the calculated result with the value of the chosen option.
    # A tolerance is used because the options are given with "approximately" (~) signs.
    tolerance = 0.01
    if math.isclose(calculated_ratio, llm_answer_value, rel_tol=tolerance):
        # The reasoning in the provided answer is also checked manually and found to be correct.
        # It correctly derives P₂/P₁ = (5/7)³ and matches it to option C.
        return "Correct"
    else:
        return (f"The final answer is incorrect. "
                f"The correct derivation is P₂/P₁ = (Δλ₁/Δλ₂)³ = (5/7)³ ≈ {calculated_ratio:.4f}. "
                f"The provided answer chose option {llm_answer_choice} with a value of {llm_answer_value}, which does not match the calculated result.")

# Execute the check
result = check_answer()
print(result)