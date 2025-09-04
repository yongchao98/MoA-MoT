import math

def check_exoplanet_period_ratio():
    """
    This function checks the correctness of the LLM's answer by recalculating the ratio of the orbital periods.
    """
    
    # Given values from the question
    delta_lambda_1 = 5  # miliangstrom
    delta_lambda_2 = 7  # miliangstrom
    
    # The LLM's chosen answer is A, which corresponds to a value of ~0.36
    llm_answer_choice = 'A'
    llm_answer_value = 0.36

    # --- Step 1: Relate wavelength shift ratio to radial velocity ratio ---
    # The amplitude of the wavelength shift (Δλ) is directly proportional to the 
    # star's radial velocity amplitude (K). So, K₁/K₂ = Δλ₁/Δλ₂.
    K1_over_K2 = delta_lambda_1 / delta_lambda_2

    # --- Step 2: Relate radial velocity ratio to semi-major axis ratio ---
    # The formula for K is K ∝ (m_p * sin(i) / sqrt(M_star)) * (1 / sqrt(a)).
    # The problem states M_star1 = M_star2 and m_p1 ≈ m_p2.
    # For the problem to be solvable, we must assume m_p1 = m_p2 and sin(i1) = sin(i2).
    # Under these assumptions, K ∝ 1/√a.
    # Therefore, K₁/K₂ = √a₂/√a₁.
    
    # --- Step 3: Calculate the ratio of the semi-major axes (a₂/a₁) ---
    # From the previous steps, Δλ₁/Δλ₂ = √a₂/√a₁.
    # So, a₂/a₁ = (Δλ₁/Δλ₂)²
    a2_over_a1 = (delta_lambda_1 / delta_lambda_2)**2
    
    # --- Step 4: Use Kepler's Third Law to find the period ratio (T₂/T₁) ---
    # Kepler's Third Law: T² ∝ a³/M_star.
    # Since M_star1 = M_star2, this simplifies to T² ∝ a³.
    # Therefore, (T₂/T₁)² = (a₂/a₁)³.
    # Taking the square root gives T₂/T₁ = (a₂/a₁)^(3/2).
    
    # --- Step 5: Final Calculation ---
    # Substitute the expression for a₂/a₁ into the period ratio equation.
    # T₂/T₁ = ( (Δλ₁/Δλ₂)² )^(3/2) = (Δλ₁/Δλ₂)³
    calculated_T2_over_T1 = (delta_lambda_1 / delta_lambda_2)**3
    
    # --- Verification ---
    # Check if the calculated result matches the value of the chosen answer option.
    # We use math.isclose to account for potential rounding in the answer choices.
    # A relative tolerance of 5% is reasonable for multiple-choice questions.
    if llm_answer_choice != 'A':
        return f"Incorrect. The LLM chose answer {llm_answer_choice}, but the correct answer is A (~0.36)."

    if math.isclose(calculated_T2_over_T1, llm_answer_value, rel_tol=0.05):
        return "Correct"
    else:
        return (f"Incorrect. The calculated ratio T₂/T₁ is {calculated_T2_over_T1:.4f}. "
                f"The value for the chosen answer 'A' is {llm_answer_value}. "
                f"While the values are close, the LLM's choice might be based on rounding. "
                f"However, the derivation logic presented by the LLM is sound and leads to this result.")

# Run the check
result = check_exoplanet_period_ratio()
print(result)