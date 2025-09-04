import math

def check_answer():
    """
    This function checks the correctness of the provided answer to the exoplanet temperature ratio problem.
    """
    # --- Problem Data ---
    # Planet 1
    m1 = 7.0  # Earth masses
    d_lambda1 = 0.03  # Angstroms

    # Planet 2
    m2 = 5.0  # Earth masses
    d_lambda2 = 0.04  # Angstroms

    # The options provided in the question
    options = {
        'A': 1.05,
        'B': 0.98,
        'C': 1.30,
        'D': 0.53
    }
    
    # The final answer provided by the LLM to be checked.
    # The prompt shows multiple answers, but the final one to be evaluated is 'A'.
    llm_answer_key = 'A'

    # --- Calculation ---
    # The ratio of equilibrium temperatures T_eq1 / T_eq2 is derived as:
    # T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ1 / Δλ2)
    
    # Calculate the mass ratio (Planet 2 / Planet 1)
    mass_ratio = m2 / m1
    
    # Calculate the Doppler shift ratio (Planet 1 / Planet 2)
    doppler_ratio = d_lambda1 / d_lambda2
    
    # Calculate the final temperature ratio
    calculated_ratio = mass_ratio * doppler_ratio

    # --- Verification ---
    # Find the option that is numerically closest to the calculated result.
    # This is done by finding the minimum absolute difference between the calculated value and the option values.
    correct_option_key = min(options, key=lambda k: abs(options[k] - calculated_ratio))

    # Check if the LLM's answer matches the calculated correct option.
    if llm_answer_key == correct_option_key:
        return "Correct"
    else:
        reason = (
            f"Incorrect. The provided answer is {llm_answer_key}, but the correct answer is {correct_option_key}.\n\n"
            f"Reasoning:\n"
            f"1. The ratio of equilibrium temperatures (T_eq1 / T_eq2) for two planets around the same star with the same albedo is sqrt(a2 / a1), where 'a' is the semi-major axis.\n"
            f"2. The semi-major axis 'a' is related to the planet mass (Mp) and the radial velocity semi-amplitude (K) by a ∝ (Mp / K)^2.\n"
            f"3. Therefore, a2 / a1 = (M_p2 / M_p1)^2 * (K1 / K2)^2.\n"
            f"4. Substituting this into the temperature ratio gives: T_eq1 / T_eq2 = sqrt((M_p2 / M_p1)^2 * (K1 / K2)^2) = (M_p2 / M_p1) * (K1 / K2).\n"
            f"5. The radial velocity semi-amplitude K is directly proportional to the Doppler shift Δλ, so K1 / K2 = Δλ1 / Δλ2.\n"
            f"6. The final formula is: T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ1 / Δλ2).\n\n"
            f"Calculation:\n"
            f"M_p2 / M_p1 = {m2} / {m1} = {mass_ratio}\n"
            f"Δλ1 / Δλ2 = {d_lambda1} / {d_lambda2} = {doppler_ratio}\n"
            f"T_eq1 / T_eq2 = {mass_ratio:.4f} * {doppler_ratio:.2f} = {calculated_ratio:.4f}\n\n"
            f"The calculated value {calculated_ratio:.4f} is approximately 0.536, which corresponds to option D (~0.53), not option {llm_answer_key} (~{options[llm_answer_key]})."
        )
        return reason

# Run the check and print the result
result = check_answer()
print(result)