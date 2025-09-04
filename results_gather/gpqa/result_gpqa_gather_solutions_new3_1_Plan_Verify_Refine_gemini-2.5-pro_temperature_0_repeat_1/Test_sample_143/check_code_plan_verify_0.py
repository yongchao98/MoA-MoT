import math

def check_correctness_of_particle_decay_distance():
    """
    Checks the correctness of the calculated mean decay distance for the meson resonance X.

    The function uses the formula:
    L = (sqrt(E^2 - (mc^2)^2) / mc^2) * (hbar*c / Gamma)
    
    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the error.
    """
    # --- Given values from the question ---
    E_GeV = 8.0         # Total energy in GeV
    m_c2_GeV = 1.2      # Rest mass energy in GeV
    Gamma_MeV = 320.0   # Decay width in MeV

    # --- Physical constants ---
    # Using a precise value for hbar*c in MeV*fm (1 fm = 1e-15 m)
    hbar_c_MeV_fm = 197.3269804

    # --- Unit Conversion ---
    # To maintain precision and consistency, convert all energy units to MeV.
    E_MeV = E_GeV * 1000
    m_c2_MeV = m_c2_GeV * 1000
    
    # --- Calculation Steps ---
    # 1. Calculate the momentum term pc = sqrt(E^2 - (mc^2)^2)
    try:
        # The term inside the square root must be non-negative.
        pc_squared_MeV2 = E_MeV**2 - m_c2_MeV**2
        if pc_squared_MeV2 < 0:
            return (f"Incorrect: The total energy E ({E_GeV} GeV) cannot be less than the rest mass energy "
                    f"m*c^2 ({m_c2_GeV} GeV). This leads to an imaginary momentum.")
        pc_MeV = math.sqrt(pc_squared_MeV2)
    except ValueError:
        # This case should be caught by the check above, but is included for robustness.
        return "Incorrect: A math error occurred during momentum calculation, likely due to E < mc^2."

    # 2. Calculate the relativistic factor gamma*beta = pc / mc^2
    gamma_beta = pc_MeV / m_c2_MeV

    # 3. Calculate the proper decay length c*tau = hbar*c / Gamma
    c_tau_fm = hbar_c_MeV_fm / Gamma_MeV

    # 4. Calculate the mean decay distance L in femtometers (fm)
    L_fm = gamma_beta * c_tau_fm

    # 5. Convert the final result to meters (1 fm = 1e-15 m)
    calculated_L_m = L_fm * 1e-15

    # --- Verification ---
    # The provided answer selects option B, which corresponds to the value 4.0655e-15 m.
    expected_L_m = 4.0655e-15

    # Check if the calculated value is close to the expected value.
    # A small tolerance (e.g., 0.5%) accounts for potential rounding differences
    # in the physical constants used to generate the problem's options.
    tolerance = 0.005  # 0.5%
    
    if abs(calculated_L_m - expected_L_m) / expected_L_m < tolerance:
        return "Correct"
    else:
        return (f"Incorrect: The calculated mean decay distance is {calculated_L_m:.4e} m, "
                f"while the value from the selected answer is {expected_L_m:.4e} m. "
                f"The relative difference is {abs(calculated_L_m - expected_L_m) / expected_L_m:.2%}, "
                f"which is outside the tolerance of {tolerance:.1%}. The calculation is sound, but the "
                f"result does not match the chosen option's value precisely, though it is very close.")

# The code block to be executed:
print(check_correctness_of_particle_decay_distance())