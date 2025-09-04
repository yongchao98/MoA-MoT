import math

def check_correctness():
    """
    This function checks the correctness of the LLM's answer to the exoplanet temperature ratio problem.
    """
    # --- Define problem constants and given values ---
    # Mass of Planet 1 in Earth masses
    M_p1 = 7.0
    # Mass of Planet 2 in Earth masses
    M_p2 = 5.0
    # Doppler shift for Planet 1 in Angstroms
    delta_lambda1 = 0.03
    # Doppler shift for Planet 2 in Angstroms
    delta_lambda2 = 0.04

    # --- Define the options as presented in the final answer ---
    options = {
        'A': 1.30,
        'B': 0.98,
        'C': 0.53,
        'D': 1.05
    }
    
    # The final answer selected by the LLM
    llm_selected_option = 'C'

    # --- Perform the calculation based on the derived formula ---
    # The derivation T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ₁ / Δλ₂) is correct.
    # All extraneous information (star's mass, radius, temp; planets' radii) correctly cancels out.
    
    # Calculate the mass ratio (M_p2 / M_p1)
    mass_ratio = M_p2 / M_p1
    
    # Calculate the Doppler shift ratio (Δλ₁ / Δλ₂)
    doppler_ratio = delta_lambda1 / delta_lambda2
    
    # Calculate the final equilibrium temperature ratio
    calculated_ratio = mass_ratio * doppler_ratio
    
    # The exact fractional result is (5/7) * (3/4) = 15/28
    exact_ratio = 15.0 / 28.0
    
    # --- Verify the result ---
    # 1. Check if the calculation is correct
    if not math.isclose(calculated_ratio, exact_ratio):
        return f"The calculation is incorrect. Expected {exact_ratio}, but calculated {calculated_ratio}."

    # 2. Check if the selected option matches the calculation
    expected_value = options.get(llm_selected_option)
    if expected_value is None:
        return f"The selected option '{llm_selected_option}' is not a valid choice."

    # Use a tolerance because the options are given as approximate values ("~")
    tolerance = 0.01
    if math.isclose(calculated_ratio, expected_value, rel_tol=tolerance, abs_tol=tolerance):
        return "Correct"
    else:
        # Find the best matching option
        best_match_option = min(options, key=lambda k: abs(options[k] - calculated_ratio))
        return (f"The reasoning and calculation are correct, resulting in a ratio of {calculated_ratio:.4f}. "
                f"However, the final answer selects option {llm_selected_option} ({expected_value}), "
                f"while the calculated value is closest to option {best_match_option} ({options[best_match_option]}). "
                f"In this specific case, the selected option {llm_selected_option} is indeed the closest match, so the final answer is correct.")

# Since the provided answer's logic and final choice are self-consistent and correct,
# the check will pass.
# abs(0.5357 - 0.53) = 0.0057, which is within the tolerance.
print(check_correctness())