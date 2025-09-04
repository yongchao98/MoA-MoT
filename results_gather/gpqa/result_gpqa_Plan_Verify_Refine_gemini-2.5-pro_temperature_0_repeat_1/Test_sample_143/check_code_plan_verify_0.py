import numpy as np

def check_meson_decay_distance():
    """
    This function checks the calculation of the mean decay distance for a meson resonance.
    It recalculates the value based on the given physical quantities and compares it
    to the value from the selected option C.
    """
    # --- Given values from the question ---
    E_X_GeV = 8.0      # Total energy in GeV
    m_X_GeV = 1.2      # Rest mass in GeV/c^2
    Gamma_X_MeV = 320.0  # Decay width in MeV

    # --- Physical Constants ---
    # hbar*c is a useful constant in particle physics, approximately 197.327 MeV*fm
    hbar_c_MeV_fm = 197.3269804
    # Conversion factor from fm to m
    fm_to_m = 1e-15

    # --- LLM's selected answer ---
    # Option C is 4.0655 * 10^-15 m
    llm_answer_value = 4.0655e-15

    # --- Calculation ---
    # For consistency, convert all energy units to MeV.
    E_X_MeV = E_X_GeV * 1000
    m_X_MeV = m_X_GeV * 1000

    # Check if the particle can exist (Energy must be >= mass)
    if E_X_MeV < m_X_MeV:
        return f"Incorrect. The total energy ({E_X_GeV} GeV) cannot be less than the rest mass ({m_X_GeV} GeV)."

    # 1. Calculate the momentum (p*c) using the relativistic energy-momentum relation:
    # E^2 = (p*c)^2 + (m*c^2)^2  =>  p*c = sqrt(E^2 - (m*c^2)^2)
    # The result will be in MeV.
    pc_MeV = np.sqrt(E_X_MeV**2 - m_X_MeV**2)

    # 2. Calculate the mean decay distance L using the formula:
    # L = (p*c / m*c^2) * (hbar*c / Gamma)
    # The units will be: (MeV / MeV) * (MeV*fm / MeV) = fm
    L_fm = (pc_MeV / m_X_MeV) * (hbar_c_MeV_fm / Gamma_X_MeV)

    # 3. Convert the result from femtometers (fm) to meters (m)
    L_m = L_fm * fm_to_m

    # --- Verification ---
    # Compare the calculated result with the LLM's answer using a small tolerance
    # for floating-point inaccuracies.
    if np.isclose(L_m, llm_answer_value, rtol=1e-5):
        return "Correct"
    else:
        return (f"Incorrect. The calculated mean decay distance is {L_m:.4e} m, "
                f"which does not match the answer from option C ({llm_answer_value:.4e} m).\n"
                f"Calculation details:\n"
                f"  - Momentum (pc): sqrt({E_X_MeV}^2 - {m_X_MeV}^2) = {pc_MeV:.2f} MeV\n"
                f"  - Decay length in fm: ({pc_MeV:.2f} / {m_X_MeV}) * ({hbar_c_MeV_fm} / {Gamma_X_MeV}) = {L_fm:.4f} fm\n"
                f"  - Decay length in m: {L_fm:.4f} fm * 1e-15 m/fm = {L_m:.4e} m")

# Execute the check
result = check_meson_decay_distance()
print(result)