import math

def check_answer():
    """
    This function calculates the ratio of equilibrium temperatures for the two planets
    and determines the correct option from the given choices.
    """
    # --- Input Parameters ---
    # Masses of the planets (in Earth masses, units cancel out)
    m_p1 = 7.0
    m_p2 = 5.0

    # Maximum Doppler shifts (in Angstroms, units cancel out)
    delta_lambda1 = 0.03
    delta_lambda2 = 0.04

    # Multiple choice options
    options = {
        "A": 0.53,
        "B": 1.30,
        "C": 0.98,
        "D": 1.05
    }

    # --- Physics Derivations ---
    # 1. The equilibrium temperature (T_eq) of a planet is proportional to 1/sqrt(a),
    #    where 'a' is the semi-major axis. All other factors (star's temperature,
    #    star's radius, albedo) are the same for both planets and cancel out in the ratio.
    #    Therefore, T_eq1 / T_eq2 = sqrt(a2 / a1).

    # 2. The radial velocity semi-amplitude (K) of the star is proportional to the
    #    planet's mass (M_p) and inversely proportional to the square root of the
    #    semi-major axis (a). K ∝ M_p / sqrt(a).

    # 3. The Doppler shift (Δλ) is directly proportional to the radial velocity (K).
    #    So, we can substitute K with Δλ in the proportionality: Δλ ∝ M_p / sqrt(a).

    # 4. From step 3, we can rearrange to express sqrt(a) in terms of M_p and Δλ:
    #    sqrt(a) ∝ M_p / Δλ.

    # 5. Now we can find the ratio sqrt(a2 / a1):
    #    sqrt(a2) / sqrt(a1) = (M_p2 / Δλ2) / (M_p1 / Δλ1)
    #    sqrt(a2 / a1) = (M_p2 / M_p1) * (Δλ1 / Δλ2)

    # 6. Combining step 1 and step 5 gives the final formula for the temperature ratio:
    #    T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ1 / Δλ2)

    # --- Calculation ---
    mass_ratio = m_p2 / m_p1
    doppler_ratio = delta_lambda1 / delta_lambda2
    temperature_ratio = mass_ratio * doppler_ratio

    # --- Find the best matching option ---
    # Find the option with the minimum absolute difference from the calculated value.
    best_option_key = min(options, key=lambda k: abs(options[k] - temperature_ratio))
    
    # --- Print Results ---
    print(f"Mass ratio (M_p2 / M_p1): {mass_ratio:.4f}")
    print(f"Doppler shift ratio (Δλ1 / Δλ2): {doppler_ratio:.4f}")
    print(f"Calculated temperature ratio (T_eq1 / T_eq2): {temperature_ratio:.4f}")
    print(f"The closest option is {best_option_key}: {options[best_option_key]}")

    # --- Verification Logic ---
    # The provided answers are a mix, but the correct calculation leads to ~0.536.
    # Let's check if the calculated best option is 'A'.
    if best_option_key == 'A':
        return "Correct"
    else:
        # This part of the code would be reached if the calculation resulted in a different option.
        # It also serves to explain why other answers are incorrect.
        return (f"Incorrect. The calculated ratio is {temperature_ratio:.4f}, "
                f"which corresponds to option {best_option_key} ({options[best_option_key]}), not the provided answer.")

# Execute the check
# The code will calculate the result and compare it to the expected option 'A'.
# If the calculation yields a result closest to option A, it will return "Correct".
# Otherwise, it will return an error message.
# This script validates the reasoning that leads to option A.
# Many of the provided LLM answers perform the calculation correctly but fail to select the right option letter.
# For example, Answer 1 calculates ~0.536 but chooses C. Answer 2 calculates ~0.536 but chooses D.
# This code confirms the calculation and the choice of option A.

result = check_answer()
print(f"\nVerification Result: The logic leading to option A is {result}.")
