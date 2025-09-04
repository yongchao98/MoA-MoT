import math

def check_exoplanet_temperature_ratio():
    """
    This function verifies the logic and calculation for the ratio of equilibrium temperatures
    between two exoplanets.

    The key physical relationships are:
    1. Equilibrium Temperature (T_eq) for a planet around a given star with a given albedo is
       proportional to 1 / sqrt(a), where 'a' is the semi-major axis.
       Therefore, T_eq1 / T_eq2 = sqrt(a2 / a1).

    2. Radial Velocity semi-amplitude (K) for a planet in a circular orbit is proportional to
       the planet's mass (Mp) and inversely proportional to the square root of the semi-major axis (a).
       K ∝ Mp / sqrt(a). This can be rearranged to a ∝ (Mp / K)^2.

    3. The radial velocity (K) is directly proportional to the observed Doppler shift (Δλ),
       since K/c = Δλ/λ₀. Thus, the ratio K1 / K2 is equal to Δλ1 / Δλ2.
    """

    # --- Given data from the question ---
    # Planet 1
    Mp1 = 7.0  # Earth masses
    delta_lambda1 = 0.03  # Angstroms

    # Planet 2
    Mp2 = 5.0  # Earth masses
    delta_lambda2 = 0.04  # Angstroms

    # --- LLM's Answer ---
    llm_answer_choice = 'C'
    llm_answer_value = 0.53
    llm_reasoning_value = 15.0 / 28.0

    # --- Verification ---

    # Step 1: Derive the formula for the temperature ratio.
    # T_eq1 / T_eq2 = sqrt(a2 / a1)
    # Substitute a ∝ (Mp / K)^2:
    # T_eq1 / T_eq2 = sqrt( (Mp2/K2)^2 / (Mp1/K1)^2 )
    # T_eq1 / T_eq2 = (Mp2 / K2) / (Mp1 / K1) = (Mp2 / Mp1) * (K1 / K2)
    # Substitute K1/K2 = delta_lambda1 / delta_lambda2:
    # T_eq1 / T_eq2 = (Mp2 / Mp1) * (delta_lambda1 / delta_lambda2)
    # This confirms the formula used in the LLM's reasoning is correct.

    # Step 2: Calculate the ratio using the derived formula.
    calculated_ratio = (Mp2 / Mp1) * (delta_lambda1 / delta_lambda2)

    # Step 3: Compare the calculated value with the LLM's reasoning.
    if not math.isclose(calculated_ratio, llm_reasoning_value, rel_tol=1e-9):
        return (f"The reasoning contains a calculation error. "
                f"The ratio should be (Mp2/Mp1) * (Δλ1/Δλ2) = ({Mp2}/{Mp1}) * ({delta_lambda1}/{delta_lambda2}) = {calculated_ratio:.4f}. "
                f"The LLM's value of 15/28 = {llm_reasoning_value:.4f} is correct, but this check failed, indicating a potential code issue.")

    # Step 4: Check if the chosen option is the closest to the calculated value.
    options = {'A': 0.98, 'B': 1.05, 'C': 0.53, 'D': 1.30}
    
    distances = {option: abs(value - calculated_ratio) for option, value in options.items()}
    closest_option = min(distances, key=distances.get)

    if closest_option != llm_answer_choice:
        return (f"The final answer is incorrect. The calculated ratio is {calculated_ratio:.4f}. "
                f"The closest option is '{closest_option}' ({options[closest_option]}), not '{llm_answer_choice}' ({llm_answer_value}).")

    # Step 5: Final check on the approximation. The value 0.53 is a rounded version of 15/28 (~0.5357).
    # The choice of C is appropriate as it's the best numerical approximation among the options.
    
    return "Correct"

# Run the check
result = check_exoplanet_temperature_ratio()
print(result)