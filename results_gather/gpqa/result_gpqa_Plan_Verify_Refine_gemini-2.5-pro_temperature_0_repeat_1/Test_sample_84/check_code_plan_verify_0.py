import math

def check_solution():
    """
    This function verifies the LLM's answer by recalculating the ratio of equilibrium temperatures.
    """
    # --- Define given values from the question ---
    # Planet 1
    m1 = 7.0  # Earth masses
    delta_lambda1 = 0.03  # Angstroms

    # Planet 2
    m2 = 5.0  # Earth masses
    delta_lambda2 = 0.04  # Angstroms

    # --- Define the LLM's answer ---
    llm_answer_option = 'B'
    llm_answer_value = 0.53

    # --- Physics Derivation ---
    # 1. The equilibrium temperature (T_eq) of a planet is proportional to 1/sqrt(a),
    #    where 'a' is the semi-major axis of its orbit. This assumes the star's luminosity
    #    and the planet's albedo are constant or equal, which is given.
    #    Therefore, the ratio of temperatures is: T_eq1 / T_eq2 = sqrt(a2 / a1).

    # 2. The radial velocity amplitude (K) induced on a star by a planet in a circular orbit
    #    is proportional to M_p / sqrt(a), where M_p is the planet's mass.
    #    (K ≈ (M_p / M_star^(1/2)) * sqrt(G/a)). Since G and M_star are constant for the system,
    #    K ∝ M_p / sqrt(a).

    # 3. We can rearrange the relation from step 2 to solve for sqrt(a):
    #    sqrt(a) ∝ M_p / K.
    #    So, the ratio sqrt(a2) / sqrt(a1) = (M_p2 / K2) / (M_p1 / K1) = (M_p2 / M_p1) * (K1 / K2).

    # 4. Combining step 1 and 3:
    #    T_eq1 / T_eq2 = sqrt(a2 / a1) = (M_p2 / M_p1) * (K1 / K2).

    # 5. The radial velocity amplitude K is directly proportional to the maximum Doppler shift Δλ
    #    (K = c * Δλ / λ₀). Since c and λ₀ are constants, K ∝ Δλ.
    #    Therefore, K1 / K2 = Δλ1 / Δλ2.

    # 6. Substituting step 5 into step 4 gives the final formula:
    #    T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ1 / Δλ2).

    # --- Calculation ---
    try:
        # Calculate the ratio of planet masses (M_p2 / M_p1)
        mass_ratio = m2 / m1

        # Calculate the ratio of Doppler shifts (Δλ1 / Δλ2)
        shift_ratio = delta_lambda1 / delta_lambda2

        # Calculate the final temperature ratio
        calculated_temp_ratio = mass_ratio * shift_ratio
    except ZeroDivisionError:
        return "Error: Division by zero in calculation. Check input values."
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Verification ---
    # The LLM's logic follows the same derivation, resulting in the same formula.
    # The LLM's code calculates: (5/7) * (0.03/0.04) = 15/28 ≈ 0.5357
    # This value is approximately 0.53, which corresponds to option B.

    # Check if the calculated value is close to the value of the chosen option.
    # The options are A) ~1.05, B) ~0.53, C) ~1.30, D) ~0.98.
    # The calculated value is ~0.536.
    if math.isclose(calculated_temp_ratio, llm_answer_value, rel_tol=0.05):
        # The LLM correctly identified that the star's properties and planets' radii were not needed for the ratio.
        # The derivation, calculation, and final choice are all consistent and correct.
        return "Correct"
    else:
        return (f"Incorrect. The calculated temperature ratio is {calculated_temp_ratio:.4f}. "
                f"The LLM's answer is option {llm_answer_option} ({llm_answer_value}), which does not match the calculation. "
                f"The correct answer should be approximately {calculated_temp_ratio:.2f}, which corresponds to option B.")

# Execute the check and print the result
result = check_solution()
print(result)