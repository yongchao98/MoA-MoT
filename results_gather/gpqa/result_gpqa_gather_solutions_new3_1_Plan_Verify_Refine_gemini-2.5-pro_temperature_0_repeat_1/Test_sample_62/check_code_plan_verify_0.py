import math

def check_correctness():
    """
    This function checks the correctness of the final answer by recalculating the photon's momentum
    based on the problem's parameters and physical principles.
    """
    # --- Given Parameters from the Question ---
    Mx_amu = 20.0  # mass of atom X in amu
    My_amu = 2.0   # mass of atom Y in amu
    R_angstrom = 2.0 # molecular bond length in angstroms
    omega = 4.0e14 # angular frequency of vibration in rad/s

    # --- Physical Constants (using precise values from CODATA 2018) ---
    amu_to_kg = 1.66053906660e-27 # kg
    h_bar = 1.054571817e-34      # J*s (reduced Planck constant)
    c = 299792458                # m/s (speed of light)

    # --- Unit Conversions ---
    R_m = R_angstrom * 1e-10 # convert angstroms to meters

    # --- Step 1: Calculate Reduced Mass (mu) ---
    # The formula for reduced mass is mu = (m1 * m2) / (m1 + m2)
    mu_amu = (Mx_amu * My_amu) / (Mx_amu + My_amu)
    mu_kg = mu_amu * amu_to_kg

    # --- Step 2: Calculate Moment of Inertia (I) ---
    # The formula for moment of inertia of a diatomic molecule is I = mu * R^2
    I = mu_kg * (R_m ** 2)

    # --- Step 3: Calculate Energy Difference (Delta_E) ---
    # The transition is from the ground state (v=0, J=0) to the first excited
    # rovibrational state (v=1, J=1).
    # The energy difference is Delta_E = h_bar*omega + h_bar^2 / I
    delta_E = (h_bar * omega) + (h_bar**2 / I)

    # --- Step 4: Calculate Photon Momentum (p) ---
    # The energy-momentum relation for a photon is E = pc, so p = E / c
    calculated_p = delta_E / c

    # --- Step 5: Compare with the chosen answer ---
    # The final answer provided is <<<B>>>, which corresponds to p = 1.4*10^(-28) N*s
    chosen_answer_value = 1.4e-28

    # We use a relative tolerance to account for potential rounding in the options.
    # A 5% tolerance is reasonable.
    if math.isclose(calculated_p, chosen_answer_value, rel_tol=0.05):
        return "Correct"
    else:
        return (f"Incorrect. The calculated momentum is {calculated_p:.4e} N·s. "
                f"The chosen answer 'B' corresponds to a momentum of {chosen_answer_value:.4e} N·s. "
                f"The calculated value does not match the chosen option's value within a 5% tolerance.")

# Execute the check and print the result
result = check_correctness()
print(result)