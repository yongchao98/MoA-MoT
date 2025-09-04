import math

def check_decay_distance():
    """
    Calculates the mean decay distance of a meson resonance and checks it against the proposed answer.
    """
    # --- Given values from the question ---
    E_X_GeV = 8.0  # Production energy in GeV
    m_X_GeV = 1.2  # Mass in GeV (representing rest energy m*c^2)
    Gamma_X_MeV = 320.0  # Decay width in MeV

    # --- Physical Constants ---
    # Using a high-precision value for h_bar * c in MeV*fm
    # 1 fm (femtometer) = 1e-15 m
    h_bar_c_MeV_fm = 197.3269804

    # --- Proposed Answer ---
    # The provided answer is C, which corresponds to 4.0655 * 10^-15 m
    proposed_answer_value = 4.0655e-15

    # --- Calculation ---
    # Step 1: Ensure consistent units. Convert all energy-related values to MeV.
    E_X_MeV = E_X_GeV * 1000
    m_X_MeV = m_X_GeV * 1000
    Gamma_X_MeV_val = Gamma_X_MeV

    # Step 2: Calculate the momentum term (pc) from the relativistic energy-momentum relation.
    # E^2 = (pc)^2 + (mc^2)^2  =>  pc = sqrt(E^2 - (mc^2)^2)
    try:
        pc_squared_MeV2 = E_X_MeV**2 - m_X_MeV**2
        if pc_squared_MeV2 < 0:
            return "Incorrect. The particle's total energy (8 GeV) cannot be less than its rest mass energy (1.2 GeV)."
        pc_MeV = math.sqrt(pc_squared_MeV2)
    except ValueError:
        return "Incorrect. Mathematical error during momentum calculation."

    # Step 3: Calculate the mean decay distance (L) using the formula:
    # L = (pc / mc^2) * (ħc / Γ)
    # The result will be in femtometers (fm).
    L_fm = (pc_MeV / m_X_MeV) * (h_bar_c_MeV_fm / Gamma_X_MeV_val)

    # Step 4: Convert the result from femtometers to meters.
    L_m = L_fm * 1e-15

    # --- Verification ---
    # Check if the calculated value matches the proposed answer within a small tolerance
    # to account for potential rounding differences in the constants used to create the options.
    # A relative tolerance of 0.1% is reasonable.
    relative_tolerance = 0.001
    if abs(L_m - proposed_answer_value) / proposed_answer_value < relative_tolerance:
        return "Correct"
    else:
        return (f"Incorrect. The calculated mean decay distance is {L_m:.5e} m, "
                f"which does not match the proposed answer's value of {proposed_answer_value:.5e} m. "
                f"The physics and calculation steps in the provided analysis are correct, but the final answer should be verified against this calculation.")

# Run the check
result = check_decay_distance()
print(result)