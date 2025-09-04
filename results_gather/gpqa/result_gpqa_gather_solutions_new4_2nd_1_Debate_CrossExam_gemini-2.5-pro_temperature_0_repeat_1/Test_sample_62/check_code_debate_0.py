import math

def check_answer():
    """
    Checks the correctness of the provided answer for the diatomic molecule problem.
    """
    # --- Constants ---
    # Using values from scipy.constants for high precision
    h_bar = 1.054571817e-34  # Reduced Planck constant in J·s
    c = 2.99792458e8         # Speed of light in m/s
    amu_to_kg = 1.66053906660e-27 # Conversion factor for amu to kg

    # --- Given Parameters from the Question ---
    Mx_amu = 20.0  # Mass of atom X in amu
    My_amu = 2.0   # Mass of atom Y in amu
    R_angstrom = 2.0 # Molecular bond length in angstroms
    omega_rad_s = 4.0e14 # Angular frequency of vibration in rad/s

    # --- LLM's Answer ---
    # The final answer provided is <<<B>>>.
    llm_choice = 'B'
    options = {
        'A': 1.9e-28,
        'B': 1.4e-28,
        'C': 1.1e-27,
        'D': 2.3e-27
    }
    llm_answer_value = options[llm_choice]

    # --- Step-by-step Calculation ---

    # 1. Convert units to SI
    Mx_kg = Mx_amu * amu_to_kg
    My_kg = My_amu * amu_to_kg
    R_m = R_angstrom * 1e-10

    # 2. Calculate reduced mass (μ)
    # μ = (m1 * m2) / (m1 + m2)
    mu = (Mx_kg * My_kg) / (Mx_kg + My_kg)

    # 3. Calculate moment of inertia (I)
    # I = μ * R^2
    I = mu * R_m**2

    # 4. Determine the transition and calculate energy difference (ΔE)
    # The fundamental state is (v=0, J=0).
    # The next lowest energy state accessible by photon absorption (Δv=+1, ΔJ=±1) is (v=1, J=1).
    # E_initial = (1/2)ħω
    # E_final = (3/2)ħω + 2B = (3/2)ħω + ħ²/I
    # ΔE = E_final - E_initial = ħω + ħ²/I
    delta_E = (h_bar * omega_rad_s) + (h_bar**2 / I)

    # 5. Calculate the photon's momentum (p)
    # p = ΔE / c
    calculated_p = delta_E / c

    # --- Verification ---
    # Check if the calculated momentum matches the value of the chosen option.
    # We use a relative tolerance to account for potential rounding in the option values.
    # A 5% tolerance is generous enough.
    if math.isclose(calculated_p, llm_answer_value, rel_tol=0.05):
        return "Correct"
    else:
        # Find which option the calculation actually matches
        best_match = ''
        min_diff = float('inf')
        for key, value in options.items():
            diff = abs(calculated_p - value) / value
            if diff < min_diff:
                min_diff = diff
                best_match = key
        
        reason = (
            f"The provided answer is incorrect.\n"
            f"The chosen option is {llm_choice}, which corresponds to a momentum of {llm_answer_value:.3e} N·s.\n"
            f"However, the step-by-step calculation yields a momentum of {calculated_p:.3e} N·s.\n"
            f"The calculated value {calculated_p:.3e} N·s is closest to option {best_match} ({options[best_match]:.3e} N·s)."
        )
        return reason

# Run the check
result = check_answer()
print(result)