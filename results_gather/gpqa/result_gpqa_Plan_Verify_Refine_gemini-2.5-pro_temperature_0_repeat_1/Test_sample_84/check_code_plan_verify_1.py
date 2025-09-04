import math

def check_planetary_temperature_ratio():
    """
    This function checks the correctness of the answer to the exoplanet problem.
    
    The key derivation steps are:
    1. The ratio of equilibrium temperatures T_eq1 / T_eq2 simplifies to sqrt(a2 / a1),
       where 'a' is the orbital distance, because star properties and albedo are the same.
    2. The radial velocity semi-amplitude K is proportional to the planet mass M_p and inversely
       proportional to the square root of the orbital distance (K ∝ M_p / sqrt(a)).
    3. Rearranging gives a ∝ (M_p / K)^2.
    4. Therefore, a2 / a1 = (M_p2 / M_p1)^2 * (K1 / K2)^2.
    5. The RV semi-amplitude K is directly proportional to the Doppler shift Δλ (K ∝ Δλ).
       So, K1 / K2 = Δλ1 / Δλ2.
    6. Substituting this in, a2 / a1 = (M_p2 / M_p1)^2 * (Δλ1 / Δλ2)^2.
    7. Finally, T_eq1 / T_eq2 = sqrt(a2 / a1) = (M_p2 / M_p1) * (Δλ1 / Δλ2).
    
    This final simplified formula shows that most of the provided data (star mass, radius,
    planet radii, etc.) are distractors.
    """
    
    # Given values from the problem
    # Planet 1
    m1 = 7.0  # Earth masses
    delta_lambda1 = 0.03  # Angstroms
    
    # Planet 2
    m2 = 5.0  # Earth masses
    delta_lambda2 = 0.04  # Angstroms
    
    # The provided answer choice is B, which is ~0.53
    expected_answer_value = 0.53
    
    # Calculate the ratio using the derived formula
    # T_ratio = (m2 / m1) * (delta_lambda1 / delta_lambda2)
    calculated_ratio = (m2 / m1) * (delta_lambda1 / delta_lambda2)
    
    # The exact value is 15/28
    exact_value = 15.0 / 28.0 # ~0.5357
    
    # Check if the calculated value is consistent with the provided answer
    # We allow a small tolerance for rounding in the answer choices.
    tolerance = 0.01
    
    if abs(calculated_ratio - expected_answer_value) < tolerance:
        return "Correct"
    else:
        return (f"Incorrect. The calculated ratio of equilibrium temperatures is {calculated_ratio:.4f} "
                f"(exact value 15/28). The provided answer is {expected_answer_value}, which is not "
                f"the result of the calculation. The derivation T_eq1/T_eq2 = (M_p2/M_p1) * (Δλ1/Δλ2) "
                f"yields (5/7) * (0.03/0.04) = 15/28 ≈ 0.5357.")

# Running the check function
result = check_planetary_temperature_ratio()
print(result)