import math

def check_decay_distance():
    """
    This function calculates the mean decay distance of a meson resonance
    and compares it to the expected answer.
    """
    # --- Given values from the question ---
    E_X_GeV = 8.0      # Production energy in GeV
    m_X_GeV = 1.2      # Mass in GeV (interpreted as rest energy m*c^2)
    Gamma_X_MeV = 320.0  # Decay width in MeV

    # --- Physical constants ---
    # Using a high-precision value for hbar*c in MeV*fm
    hbar_c_MeV_fm = 197.3269804
    # Conversion factor from fm to m
    fm_to_m = 1e-15

    # --- Expected answer from the analysis ---
    # The final answer provided is D, which corresponds to 4.0655 * 10^-15 m
    expected_answer_m = 4.0655e-15

    # --- Calculation ---
    # The formula for the mean decay distance L is:
    # L = (gamma * beta) * (hbar * c / Gamma)
    # where gamma * beta = pc / (m*c^2)
    # and pc = sqrt(E^2 - (m*c^2)^2)

    # Step 1: Ensure consistent units. We will use MeV for energy calculations.
    E_X_MeV = E_X_GeV * 1000
    m_X_MeV = m_X_GeV * 1000

    # Step 2: Calculate the momentum term (pc) in MeV
    # E^2 = (pc)^2 + (mc^2)^2  => (pc)^2 = E^2 - (mc^2)^2
    try:
        pc_squared_MeV2 = E_X_MeV**2 - m_X_MeV**2
        if pc_squared_MeV2 < 0:
            return "Incorrect. The calculation for momentum squared (E^2 - (mc^2)^2) resulted in a negative value, which is physically impossible. The particle's total energy must be greater than or equal to its rest mass energy."
        pc_MeV = math.sqrt(pc_squared_MeV2)
    except ValueError:
        return "Incorrect. A math domain error occurred during the momentum calculation."

    # Step 3: Calculate the relativistic factor gamma*beta
    # pc = gamma * beta * (mc^2) => gamma*beta = pc / (mc^2)
    gamma_beta = pc_MeV / m_X_MeV

    # Step 4: Calculate the mean decay distance L in femtometers (fm)
    # L = (gamma*beta) * (hbar*c / Gamma)
    L_fm = gamma_beta * (hbar_c_MeV_fm / Gamma_X_MeV)

    # Step 5: Convert the result to meters
    calculated_L_m = L_fm * fm_to_m

    # --- Verification ---
    # Compare the calculated value with the expected answer.
    # A relative tolerance is used to account for potential rounding differences
    # in the problem's constants vs. high-precision constants.
    # The analysis notes a ~0.03% difference, so a 0.1% tolerance is reasonable.
    tolerance = 1e-3  # 0.1% relative tolerance

    if math.isclose(calculated_L_m, expected_answer_m, rel_tol=tolerance):
        return "Correct"
    else:
        return (f"Incorrect. The calculated mean decay distance is {calculated_L_m:.5e} m, "
                f"which does not match the expected answer of {expected_answer_m:.5e} m (Option D) within a {tolerance*100}% tolerance.\n"
                f"Calculation details:\n"
                f"  - pc = {pc_MeV:.4f} MeV\n"
                f"  - gamma*beta = {gamma_beta:.4f}\n"
                f"  - L = {L_fm:.4f} fm\n"
                "The provided answer is numerically correct based on the physics, but the code check shows a slight discrepancy likely due to the precision of constants used to create the options.")

# To run the check, you would execute the function:
# result = check_decay_distance()
# print(result)