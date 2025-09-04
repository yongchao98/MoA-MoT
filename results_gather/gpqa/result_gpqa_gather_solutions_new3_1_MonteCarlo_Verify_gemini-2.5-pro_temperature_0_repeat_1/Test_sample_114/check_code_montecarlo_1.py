import math

def check_answer():
    """
    Checks the correctness of the final answer for the electron-positron annihilation problem.
    """
    # Given Lorentz factors
    gamma_e = 4
    gamma_p = 2

    # The options provided in the question
    options = {'A': 96, 'B': 172, 'C': 74, 'D': 138}
    
    # The final answer provided by the LLM
    llm_answer_key = 'D'

    # --- Step 1: Derive the formula for the angle ---
    # This is based on conservation of energy and momentum.
    # Let m=1, c=1 for simplicity (as they cancel out).
    # E_initial = (gamma_e + gamma_p)
    # p_initial = sqrt(gamma_e^2 - 1) - sqrt(gamma_p^2 - 1)
    # E_final = 2 * E_photon => E_photon = (gamma_e + gamma_p) / 2
    # p_final = 2 * p_photon * cos(theta) = 2 * E_photon * cos(theta)
    # p_initial = p_final => sqrt(gamma_e^2 - 1) - sqrt(gamma_p^2 - 1) = (gamma_e + gamma_p) * cos(theta)
    # cos(theta) = (sqrt(gamma_e^2 - 1) - sqrt(gamma_p^2 - 1)) / (gamma_e + gamma_p)

    # --- Step 2: Calculate the angle ---
    try:
        cos_theta = (math.sqrt(gamma_e**2 - 1) - math.sqrt(gamma_p**2 - 1)) / (gamma_e + gamma_p)
        
        # Angle theta in radians
        theta_rad = math.acos(cos_theta)
        
        # The total angle between the photons is 2*theta. Convert to degrees.
        calculated_angle = 2 * math.degrees(theta_rad)
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Step 3: Check the correctness of the LLM's answer ---
    llm_answer_value = options.get(llm_answer_key)
    if llm_answer_value is None:
        return f"The provided answer key '{llm_answer_key}' is not a valid option."

    # Find which option is closest to the calculated value
    closest_option_key = min(options, key=lambda k: abs(options[k] - calculated_angle))
    
    if closest_option_key == llm_answer_key:
        return "Correct"
    else:
        return (f"Incorrect. The calculated angle is approximately {calculated_angle:.2f} degrees. "
                f"This is closest to option {closest_option_key} ({options[closest_option_key]} degrees), "
                f"but the provided answer was option {llm_answer_key} ({llm_answer_value} degrees).")

# Run the check
result = check_answer()
print(result)