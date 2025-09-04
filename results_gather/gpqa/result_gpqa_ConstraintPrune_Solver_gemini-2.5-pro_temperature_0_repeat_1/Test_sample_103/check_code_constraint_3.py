import math

def check_correctness_of_answer():
    """
    This function checks the correctness of the provided answer to the exoplanet problem.

    Problem Summary:
    - Two identical stars (M_star1 = M_star2).
    - Two identical planets (m_p1 = m_p2).
    - Circular orbits.
    - Wavelength shift for system 1 (Δλ₁) = 5 mÅ.
    - Wavelength shift for system 2 (Δλ₂) = 7 mÅ.
    - Question: What is the ratio of orbital periods, T₂ / T₁?
    - Provided Answer: C (~0.36)

    Physics Derivations:
    1.  The observed wavelength shift (Δλ) is directly proportional to the star's radial velocity amplitude (K) from the Doppler effect.
        Therefore, K ∝ Δλ, which implies K₂ / K₁ = Δλ₂ / Δλ₁.

    2.  From the conservation of momentum in a two-body system (M_star * v_star = m_p * v_p), and noting that the star's velocity amplitude K is its orbital speed v_star, we get K = (m_p / M_star) * v_p.
        Since M_star and m_p are the same for both systems, K is directly proportional to the planet's orbital velocity (v_p).
        Therefore, K ∝ v_p, which implies v_p₂ / v_p₁ = K₂ / K₁.

    3.  Combining the first two points: v_p₂ / v_p₁ = Δλ₂ / Δλ₁.

    4.  For a circular orbit, the planet's orbital period (T) is related to its orbital velocity (v_p).
        - From orbital mechanics, v_p = sqrt(G * M_star / a), where 'a' is the orbital radius.
        - The period is T = 2πa / v_p.
        - From these, we can derive that T is inversely proportional to the cube of the planet's velocity (T ∝ 1 / v_p³), given that the star's mass (M_star) is constant.

    5.  Therefore, the ratio of the periods is:
        T₂ / T₁ = (1 / v_p₂³) / (1 / v_p₁³) = (v_p₁ / v_p₂)³ = (Δλ₁ / Δλ₂)³.
    """

    # Given values from the problem
    delta_lambda_1 = 5.0
    delta_lambda_2 = 7.0

    # The answer to check is 'C', which corresponds to the value ~0.36
    llm_answer_choice = 'C'
    options = {'A': 0.85, 'B': 1.96, 'C': 0.36, 'D': 1.40}
    
    if llm_answer_choice not in options:
        return f"The provided answer '{llm_answer_choice}' is not a valid option."

    llm_answer_value = options[llm_answer_choice]

    # Calculate the expected ratio of the periods based on the physics
    # T₂ / T₁ = (Δλ₁ / Δλ₂)³
    expected_ratio = (delta_lambda_1 / delta_lambda_2) ** 3

    # Check if the value from the provided answer is close to the calculated value.
    # A relative tolerance of 5% is reasonable for approximate multiple-choice answers.
    if math.isclose(expected_ratio, llm_answer_value, rel_tol=0.05):
        return "Correct"
    else:
        # Find the closest option to the calculated result
        closest_option = min(options.keys(), key=lambda k: abs(options[k] - expected_ratio))
        
        reason = (
            f"The answer is incorrect. The provided answer is '{llm_answer_choice}' ({llm_answer_value}), "
            f"but the calculated period ratio T₂/T₁ is (5/7)³ ≈ {expected_ratio:.4f}. "
            f"This calculated value is closest to option '{closest_option}' ({options[closest_option]})."
        )
        return reason

# Execute the check and print the result
result = check_correctness_of_answer()
print(result)