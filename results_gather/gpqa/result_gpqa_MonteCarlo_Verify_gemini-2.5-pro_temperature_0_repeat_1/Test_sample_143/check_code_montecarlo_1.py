import math

def check_decay_distance():
    """
    Calculates the mean decay distance of a meson resonance and checks it against the provided answer.
    """
    # --- Constants ---
    # hbar * c, a useful constant in particle physics. Value from CODATA 2018.
    # Units: MeV * fm (Mega-electron-volt * femtometer)
    # 1 fm = 1e-15 m
    hbar_c_MeV_fm = 197.3269804

    # --- Given parameters from the question ---
    E_X_GeV = 8.0      # Production energy in GeV
    m_X_GeV = 1.2      # Mass in GeV/c^2
    Gamma_X_MeV = 320.0  # Decay width in MeV

    # --- Convert all units to be consistent (MeV) ---
    E_X_MeV = E_X_GeV * 1000
    m_X_c2_MeV = m_X_GeV * 1000

    # --- Constraint Check ---
    # The total energy must be greater than the rest mass energy for the particle to exist and move.
    if E_X_MeV <= m_X_c2_MeV:
        return (f"Incorrect: The production energy E_X ({E_X_MeV} MeV) must be greater than "
                f"the rest mass energy m_X*c^2 ({m_X_c2_MeV} MeV).")

    # --- Calculation ---
    # The formula for the mean decay distance (d) is: d = (p*c / m*c^2) * (hbar*c / Gamma)
    # where p is the momentum.

    # 1. Calculate p*c (momentum * speed of light) using the relativistic energy-momentum relation:
    # E^2 = (p*c)^2 + (m*c^2)^2  =>  p*c = sqrt(E^2 - (m*c^2)^2)
    pc_MeV = math.sqrt(E_X_MeV**2 - m_X_c2_MeV**2)

    # 2. Calculate the mean decay distance in femtometers (fm)
    mean_decay_distance_fm = (pc_MeV / m_X_c2_MeV) * (hbar_c_MeV_fm / Gamma_X_MeV)

    # 3. Convert the final distance to meters (1 fm = 1e-15 m)
    calculated_distance_m = mean_decay_distance_fm * 1e-15

    # --- Answer Verification ---
    # The provided answer is option A.
    llm_answer_option = 'A'
    options = {
        'A': 4.0655e-15,
        'B': 5.0223e-16,
        'C': 5.0223e-15,
        'D': 4.0655e-16
    }
    
    # Get the numerical value of the given answer
    llm_answer_value = options[llm_answer_option]

    # --- Comparison ---
    # Check if the calculated value is close to the answer's value.
    # A relative tolerance of 1% is reasonable to account for different rounding of constants.
    tolerance = 0.01 
    
    if math.isclose(calculated_distance_m, llm_answer_value, rel_tol=tolerance):
        return "Correct"
    else:
        relative_error = abs(calculated_distance_m - llm_answer_value) / llm_answer_value
        return (f"Incorrect: The calculated mean decay distance is {calculated_distance_m:.4e} m. "
                f"The value for option {llm_answer_option} is {llm_answer_value:.4e} m. "
                f"The relative error between the calculated value and the provided answer is {relative_error:.4%}, "
                f"which is outside the acceptable tolerance of {tolerance:.2%}.")

# Execute the check
result = check_decay_distance()
print(result)