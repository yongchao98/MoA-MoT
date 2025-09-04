import numpy as np

def check_correctness():
    """
    Checks the correctness of the calculated mean decay distance for the meson resonance X.
    """
    # --- Given values from the question ---
    E_X_GeV = 8.0      # Total energy in GeV
    m_X_c2_GeV = 1.2   # Rest mass energy in GeV
    Gamma_X_MeV = 320.0 # Decay width in MeV

    # --- Physical constants ---
    # hbar*c in MeV*fm (1 fm = 1e-15 m). Using a high-precision value.
    hbar_c_MeV_fm = 197.3269804

    # --- The final answer provided by the LLM ---
    # The LLM chose option B, which corresponds to 4.0655 * 10^-15 m
    llm_answer_value = 4.0655e-15

    # --- Calculation ---
    # The formula for mean decay distance L is:
    # L = (gamma * beta) * (hbar * c / Gamma)
    # where gamma * beta = pc / mc^2 = sqrt(E^2 - (mc^2)^2) / mc^2

    # Step 1: Check if the input values are physically valid.
    if E_X_GeV <= m_X_c2_GeV:
        return "Constraint not satisfied: The total energy (E_X) must be greater than the rest mass energy (m_X_c2) for the particle to be moving."

    # Step 2: Calculate the relativistic factor (gamma * beta).
    # The units of GeV cancel out, so no conversion is needed for this part.
    gamma_beta = np.sqrt(E_X_GeV**2 - m_X_c2_GeV**2) / m_X_c2_GeV

    # Step 3: Calculate the mean decay distance in femtometers (fm).
    # Here we use Gamma in MeV and hbar*c in MeV*fm.
    L_fm = gamma_beta * (hbar_c_MeV_fm / Gamma_X_MeV)

    # Step 4: Convert the result to meters.
    calculated_L_m = L_fm * 1e-15

    # --- Verification ---
    # Check if the LLM's answer value is close to the calculated value.
    # A relative tolerance of 0.1% (1e-3) is appropriate given the precision of the options.
    if np.isclose(calculated_L_m, llm_answer_value, rtol=1e-3):
        return "Correct"
    else:
        return (f"Incorrect. The calculated mean decay distance is approximately {calculated_L_m:.4e} m. "
                f"The provided answer was {llm_answer_value:.4e} m. The values do not match.")

# Run the check
result = check_correctness()
print(result)