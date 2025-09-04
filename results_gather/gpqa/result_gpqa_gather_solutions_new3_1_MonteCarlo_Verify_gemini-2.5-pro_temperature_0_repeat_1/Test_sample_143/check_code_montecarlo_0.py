import math

def check_answer():
    """
    This function checks the correctness of the calculated mean decay distance for the meson resonance X.
    """
    # Given values from the question
    E_X = 8.0  # Production energy in GeV
    m_X = 1.2  # Mass (mc^2) in GeV
    Gamma_X = 320.0  # Decay width in MeV

    # Physical constants
    # hbar * c is a useful constant in particle physics.
    # hbar_c ≈ 197.327 MeV·fm (where 1 fm = 10^-15 m)
    hbar_c_MeV_fm = 197.327

    # The final answer from the LLM is C, which corresponds to 4.0655 * 10^-15 m
    expected_value = 4.0655e-15

    # --- Calculation ---
    # The formula for the mean decay distance (L) in the lab frame is:
    # L = βγ * cτ
    # where τ = ħ/Γ is the proper lifetime.
    # This can be rewritten as L = βγ * (ħc / Γ)

    # The term βγ can be calculated from the energy-momentum relation:
    # E² = (pc)² + (mc²)²
    # pc = βγ * mc²
    # So, βγ = pc / mc² = sqrt(E² - (mc²)²) / mc²

    # Step 1: Calculate βγ
    # Ensure all units are consistent. E_X and m_X are already in GeV.
    try:
        if E_X < m_X:
            return "Incorrect. The total energy (8 GeV) cannot be less than the rest mass energy (1.2 GeV)."
        
        # Calculate the momentum term (pc) in GeV
        pc = math.sqrt(E_X**2 - m_X**2)
        
        # Calculate βγ (dimensionless)
        beta_gamma = pc / m_X
    except ValueError:
        return "Incorrect. Calculation error: Cannot take the square root of a negative number. Energy must be greater than or equal to mass."

    # Step 2: Calculate the proper decay length (cτ)
    # cτ = ħc / Γ
    # We use ħc in MeV·fm and Γ in MeV. The MeV units will cancel.
    c_tau_fm = hbar_c_MeV_fm / Gamma_X

    # Step 3: Calculate the mean decay distance (L) in femtometers
    L_fm = beta_gamma * c_tau_fm

    # Step 4: Convert the result to meters
    # 1 fm = 10^-15 m
    calculated_L_meters = L_fm * 1e-15

    # Step 5: Compare the calculated value with the expected value from option C
    # Use a relative tolerance to account for rounding differences in constants
    if math.isclose(calculated_L_meters, expected_value, rel_tol=1e-4):
        return "Correct"
    else:
        return (f"Incorrect. The calculated mean decay distance is {calculated_L_meters:.4e} m, "
                f"which does not match the expected value of {expected_value:.4e} m from option C. "
                f"The calculation steps are: βγ = sqrt(E² - m²) / m ≈ {beta_gamma:.4f}; "
                f"cτ = ħc / Γ ≈ {c_tau_fm:.4f} fm; "
                f"L = βγ * cτ ≈ {L_fm:.4f} fm.")

# Run the check and print the result
print(check_answer())