import math

def check_answer():
    """
    Checks the correctness of the LLM's answer for the diatomic molecule problem.
    """
    # --- Define Physical Constants ---
    # Reduced Planck constant (J*s)
    hbar = 1.054571817e-34
    # Atomic Mass Unit to kg conversion
    amu_to_kg = 1.66053906660e-27
    # Speed of light (m/s)
    c = 299792458

    # --- Given Parameters from the Question ---
    # Mass of atom X (amu)
    Mx_amu = 20.0
    # Mass of atom Y (amu)
    My_amu = 2.0
    # Molecular bond length (Angstroms)
    R_angstrom = 2.0
    # Angular frequency of vibration (rad/s)
    omega = 4.0e14

    # --- Convert units to SI ---
    R_m = R_angstrom * 1e-10  # Convert Angstroms to meters

    # --- Step-by-step Calculation ---

    # 1. Calculate the reduced mass (mu) in kg
    # mu = (Mx * My) / (Mx + My)
    mu_amu = (Mx_amu * My_amu) / (Mx_amu + My_amu)
    mu_kg = mu_amu * amu_to_kg

    # 2. Calculate the moment of inertia (I)
    # I = mu * R^2
    I = mu_kg * (R_m ** 2)

    # 3. Determine the energy of the transition (delta_E)
    # The transition is from the ground state (v=0, J=0) to the first
    # accessible excited state (v=1, J=1).
    # E(v, J) = hbar*omega*(v + 1/2) + (hbar^2 / (2*I)) * J*(J+1)
    # E_initial = E(0, 0) = (1/2)*hbar*omega
    # E_final = E(1, 1) = (3/2)*hbar*omega + (hbar^2 / I)
    # delta_E = E_final - E_initial = hbar*omega + hbar^2 / I
    
    vibrational_energy_term = hbar * omega
    rotational_energy_term = (hbar ** 2) / I
    delta_E = vibrational_energy_term + rotational_energy_term

    # 4. Calculate the photon's momentum (p)
    # p = delta_E / c
    calculated_p = delta_E / c

    # --- Check the LLM's Answer ---
    # The LLM's final answer is <<<D>>>, which corresponds to p = 1.4*10^(-28) N*s
    llm_answer_value = 1.4e-28

    # Compare the calculated value with the LLM's answer value.
    # We use math.isclose() for robust floating-point comparison.
    # A relative tolerance of 5% is reasonable given potential rounding in the options.
    if math.isclose(calculated_p, llm_answer_value, rel_tol=0.05):
        return "Correct"
    else:
        # Provide a detailed reason for the error
        reason = (
            f"The answer is incorrect.\n"
            f"The final answer provided was <<<D>>>, which corresponds to a momentum of {llm_answer_value:.3e} N*s.\n"
            f"The calculation based on the problem's parameters yields a different result.\n"
            f"Step-by-step verification:\n"
            f"1. Reduced Mass (μ): {mu_kg:.4e} kg\n"
            f"2. Moment of Inertia (I): {I:.4e} kg*m^2\n"
            f"3. Transition Energy (ΔE = ħω + ħ²/I): {delta_E:.4e} J\n"
            f"4. Calculated Momentum (p = ΔE/c): {calculated_p:.4e} N*s\n"
            f"The calculated momentum {calculated_p:.4e} N*s does not match the answer's value of {llm_answer_value:.3e} N*s within a 5% tolerance."
        )
        return reason

# Run the check
result = check_answer()
print(result)