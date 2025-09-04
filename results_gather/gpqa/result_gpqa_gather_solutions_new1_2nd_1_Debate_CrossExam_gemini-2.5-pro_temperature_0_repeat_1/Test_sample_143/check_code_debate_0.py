import math

def check_correctness():
    """
    This function calculates the mean decay distance of a meson resonance X
    and compares it with the proposed answer.
    """
    # --- Define Given Values and Constants ---
    # Given values from the question
    E_X_GeV = 8.0      # Total energy in GeV
    m_X_c2_GeV = 1.2   # Rest mass energy in GeV
    Gamma_X_MeV = 320.0  # Decay width in MeV

    # Physical constant hbar*c in MeV*fm
    # Using a standard precise value: 197.3269804 MeV·fm
    hbar_c_MeV_fm = 197.3269804

    # The proposed answer from the LLM is B, which corresponds to 4.0655 * 10^-15 m
    proposed_answer_value = 4.0655e-15

    # --- Perform Calculation ---
    # Step 1: Ensure consistent units. Convert everything to GeV.
    Gamma_X_GeV = Gamma_X_MeV / 1000.0
    hbar_c_GeV_fm = hbar_c_MeV_fm / 1000.0 # Now in GeV*fm

    # Step 2: Calculate the momentum term (pc) in GeV using the relativistic energy-momentum relation.
    # E^2 = (pc)^2 + (mc^2)^2  =>  pc = sqrt(E^2 - (mc^2)^2)
    try:
        pc_squared_GeV2 = E_X_GeV**2 - m_X_c2_GeV**2
        if pc_squared_GeV2 < 0:
            return "Incorrect: The total energy (8 GeV) is less than the rest mass energy, which is physically impossible. Cannot calculate momentum."
        pc_GeV = math.sqrt(pc_squared_GeV2)
    except ValueError:
        return "Incorrect: A calculation error occurred while finding the momentum."

    # Step 3: Calculate the mean decay distance (L) using the full formula.
    # L = (pc / mc^2) * (hbar*c / Gamma)
    # The result will be in femtometers (fm) because hbar_c is in GeV*fm
    # and all other energy units (GeV) cancel out.
    if Gamma_X_GeV == 0:
        return "Incorrect: Decay width cannot be zero."
    
    L_fm = (pc_GeV / m_X_c2_GeV) * (hbar_c_GeV_fm / Gamma_X_GeV)

    # Step 4: Convert the final result from femtometers (fm) to meters (m).
    # 1 fm = 1e-15 m
    calculated_L_m = L_fm * 1e-15

    # --- Verify the Answer ---
    # Compare the calculated value with the proposed answer using a relative tolerance
    # to account for floating-point inaccuracies and different precision of constants.
    # A tolerance of 0.1% is reasonable for this type of physics problem.
    if math.isclose(calculated_L_m, proposed_answer_value, rel_tol=0.001):
        return "Correct"
    else:
        return (f"Incorrect: The calculated mean decay distance is {calculated_L_m:.4e} m. "
                f"The proposed answer is {proposed_answer_value:.4e} m. "
                f"The values do not match within a 0.1% tolerance. The calculation followed the formula L = (pc / mc^2) * (hbar*c / Gamma).")

# Run the check
result = check_correctness()
print(result)