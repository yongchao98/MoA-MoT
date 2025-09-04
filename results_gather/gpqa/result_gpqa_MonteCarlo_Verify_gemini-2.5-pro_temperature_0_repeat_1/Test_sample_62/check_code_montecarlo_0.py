import math

def check_answer():
    """
    This function checks the correctness of the provided answer for the photon momentum problem.
    It calculates the required momentum from the given physical parameters and compares it
    to the value in the selected option D.
    """
    # --- Define Constants ---
    # Using standard high-precision values for physical constants
    H_BAR = 1.054571817e-34  # Reduced Planck constant in J·s
    C = 2.99792458e8         # Speed of light in m/s
    AMU_TO_KG = 1.66053906660e-27 # Conversion factor from atomic mass units to kg

    # --- Input Parameters from the Question ---
    Mx_amu = 20.0  # Mass of atom X in amu
    My_amu = 2.0   # Mass of atom Y in amu
    R_m = 2.0e-10  # Molecular bond length in meters (2 angstroms)
    omega_rad_s = 4.0e14 # Angular frequency of vibration in rad/s

    # The answer to check is Option D
    option_d_value = 1.4e-28

    # --- Step 1: Calculate the reduced mass (μ) ---
    # Formula: μ = (Mx * My) / (Mx + My)
    mu_amu = (Mx_amu * My_amu) / (Mx_amu + My_amu)
    mu_kg = mu_amu * AMU_TO_KG

    # --- Step 2: Calculate the moment of inertia (I) ---
    # Formula: I = μ * R^2
    moment_of_inertia = mu_kg * (R_m ** 2)

    # --- Step 3: Determine the energy of the transition (ΔE) ---
    # The problem asks for the transition from the fundamental state to the next lowest energy state.
    # Fundamental state: v=0, J=0. Energy E_initial = (1/2)ħω.
    # Selection rules for photon absorption: Δv = +1, ΔJ = ±1.
    # Starting from J=0, the only possibility is ΔJ = +1.
    # So, the final state is v=1, J=1.
    # Energy of final state: E_final = (3/2)ħω + (ħ^2 / 2I) * 1*(1+1) = (3/2)ħω + ħ^2/I.
    # Energy of transition: ΔE = E_final - E_initial = (ħω + ħ^2/I).
    
    vibrational_energy_term = H_BAR * omega_rad_s
    rotational_energy_term = (H_BAR ** 2) / moment_of_inertia
    delta_E = vibrational_energy_term + rotational_energy_term

    # --- Step 4: Calculate the photon momentum (p) ---
    # Formula: p = ΔE / c
    calculated_momentum = delta_E / C

    # --- Step 5: Check the correctness of the answer ---
    # We compare the calculated momentum with the value from option D.
    # We use a relative tolerance to account for potential rounding in the options.
    # A 5% relative tolerance is reasonable for this type of problem.
    if math.isclose(calculated_momentum, option_d_value, rel_tol=0.05):
        return "Correct"
    else:
        # If the answer is incorrect, provide the calculated value for debugging.
        reason = (
            f"The answer is incorrect.\n"
            f"The calculated momentum is p = {calculated_momentum:.4e} N·s.\n"
            f"The value from option D is p = {option_d_value:.4e} N·s.\n"
            f"These values do not match within a 5% relative tolerance.\n"
            f"Intermediate values:\n"
            f"  - Reduced Mass (μ): {mu_kg:.4e} kg\n"
            f"  - Moment of Inertia (I): {moment_of_inertia:.4e} kg·m²\n"
            f"  - Transition Energy (ΔE): {delta_E:.4e} J"
        )
        return reason

# Run the check and print the result
result = check_answer()
print(result)