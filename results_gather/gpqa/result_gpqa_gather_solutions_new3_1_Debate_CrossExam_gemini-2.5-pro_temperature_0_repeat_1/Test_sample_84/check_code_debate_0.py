import math

def check_answer():
    """
    This function checks the correctness of the provided answer by recalculating the result
    based on the physical principles outlined in the problem.
    """
    # --- Given values from the question ---
    # Planet 1
    m_p1 = 7.0  # Earth masses
    delta_lambda1 = 0.03  # Angstroms

    # Planet 2
    m_p2 = 5.0  # Earth masses
    delta_lambda2 = 0.04  # Angstroms

    # --- LLM's final answer ---
    # The final answer from the LLM is C, which corresponds to ~0.53
    llm_option = 'C'
    options = {'A': 1.05, 'B': 1.30, 'C': 0.53, 'D': 0.98}
    llm_value = options[llm_option]

    # --- Calculation ---
    # The ratio of equilibrium temperatures T_eq1 / T_eq2 can be derived as:
    # T_eq1 / T_eq2 = (M_p2 / M_p1) * (K1 / K2)
    # where K is the radial velocity semi-amplitude.
    # Since K is directly proportional to the Doppler shift Δλ, we have:
    # K1 / K2 = Δλ1 / Δλ2
    # Therefore, T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ1 / Δλ2)

    # Calculate the ratio of masses (Planet 2 / Planet 1)
    mass_ratio = m_p2 / m_p1

    # Calculate the ratio of Doppler shifts (Planet 1 / Planet 2)
    shift_ratio = delta_lambda1 / delta_lambda2

    # Calculate the final temperature ratio
    temp_ratio = mass_ratio * shift_ratio
    
    # --- Verification ---
    # Check if the calculated value is close to the value of the chosen option.
    # We use a tolerance of 5% to account for rounding in the options.
    if not math.isclose(temp_ratio, llm_value, rel_tol=0.05):
        return (f"Incorrect. The calculated ratio is {temp_ratio:.4f} ({m_p2}/{m_p1} * {delta_lambda1}/{delta_lambda2} = {15/28}). "
                f"This corresponds to option C (~0.53). The provided answer was {llm_option} ({llm_value}).")

    # Check the derivation logic. Most answers use a correct derivation.
    # The key steps are:
    # 1. T_eq1 / T_eq2 = sqrt(a2 / a1)
    # 2. a ∝ (Mp/K)^2
    # 3. K ∝ Δλ
    # Combining these gives T_eq1 / T_eq2 = (Mp2/Mp1) * (Δλ1/Δλ2)
    # The provided solution follows this logic.

    return "Correct"

# Print the result of the check
print(check_answer())