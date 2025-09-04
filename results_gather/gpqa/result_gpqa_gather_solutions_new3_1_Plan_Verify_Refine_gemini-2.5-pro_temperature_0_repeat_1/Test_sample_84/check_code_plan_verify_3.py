import math

def check_planet_temperature_ratio():
    """
    Calculates the ratio of equilibrium temperatures for two planets
    based on their masses and the radial velocity shifts they induce on their star.
    """
    # --- Given values from the question ---
    # Mass of Planet 1 in Earth masses
    m_p1 = 7.0
    # Mass of Planet 2 in Earth masses
    m_p2 = 5.0
    # Doppler shift for Planet 1 in Angstroms
    delta_lambda1 = 0.03
    # Doppler shift for Planet 2 in Angstroms
    delta_lambda2 = 0.04

    # --- Derivation and Calculation ---
    # The ratio of equilibrium temperatures T_eq1 / T_eq2 is derived as:
    # T_eq1 / T_eq2 = (M_p2 / M_p1) * (K1 / K2)
    # where K is the radial velocity semi-amplitude.
    # Since K is proportional to the Doppler shift Δλ, we have:
    # K1 / K2 = Δλ₁ / Δλ₂
    # So, the final formula is:
    # T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ₁ / Δλ₂)

    # Calculate the ratio
    calculated_ratio = (m_p2 / m_p1) * (delta_lambda1 / delta_lambda2)

    # --- Options from the question ---
    options = {
        "A": 0.53,
        "B": 0.98,
        "C": 1.30,
        "D": 1.05
    }
    
    # --- Verification ---
    # Find the closest option to the calculated result
    closest_option = min(options, key=lambda k: abs(options[k] - calculated_ratio))

    # The expected answer is A, which is ~0.53
    expected_answer_letter = "A"

    print(f"Calculation Steps:")
    print(f"Mass ratio (M_p2 / M_p1) = {m_p2} / {m_p1} = {m_p2 / m_p1:.4f}")
    print(f"Doppler shift ratio (Δλ₁ / Δλ₂) = {delta_lambda1} / {delta_lambda2} = {delta_lambda1 / delta_lambda2:.4f}")
    print(f"Calculated Temperature Ratio (T_eq1 / T_eq2) = ({m_p2 / m_p1:.4f}) * ({delta_lambda1 / delta_lambda2:.4f}) = {calculated_ratio:.4f}")
    print("-" * 30)
    print(f"The calculated ratio is {calculated_ratio:.4f}.")
    print(f"The closest option is {closest_option} ({options[closest_option]}).")

    # Check if the closest option matches the expected correct answer
    if closest_option == expected_answer_letter:
        # Check if the calculated value is reasonably close to the option value
        if math.isclose(calculated_ratio, options[expected_answer_letter], rel_tol=0.05):
            return "Correct"
        else:
            return (f"Incorrect. The logic leads to option {closest_option}, but the numerical value "
                    f"{calculated_ratio:.4f} is not close enough to the option value {options[closest_option]}.")
    else:
        return (f"Incorrect. The correct answer is {closest_option} with a value of ~{calculated_ratio:.4f}. "
                f"Many of the provided LLM answers chose other options, often due to misinterpreting the final result.")

# Run the verification
print(check_planet_temperature_ratio())
