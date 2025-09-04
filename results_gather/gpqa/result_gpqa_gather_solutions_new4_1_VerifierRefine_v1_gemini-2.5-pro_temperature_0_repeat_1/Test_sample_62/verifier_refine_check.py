import math

def check_diatomic_molecule_momentum():
    """
    Calculates the momentum of a photon for a rovibrational transition
    and checks it against the provided answer.
    """
    # --- 1. Define Constants and Given Parameters ---
    # Physical constants in SI units
    h_bar = 1.054571817e-34  # Reduced Planck constant (J·s)
    c = 299792458          # Speed of light (m/s)
    amu_to_kg = 1.660539e-27 # Atomic mass unit to kg conversion

    # Parameters from the question
    Mx_amu = 20.0
    My_amu = 2.0
    R_angstrom = 2.0
    omega = 4.0e14  # rad/s

    # Convert parameters to SI units
    R_m = R_angstrom * 1e-10

    # --- 2. Calculate Reduced Mass (μ) ---
    try:
        mu_amu = (Mx_amu * My_amu) / (Mx_amu + My_amu)
        mu_kg = mu_amu * amu_to_kg
    except ZeroDivisionError:
        return "Error: Sum of masses is zero, cannot calculate reduced mass."

    # --- 3. Calculate Moment of Inertia (I) ---
    I = mu_kg * (R_m ** 2)

    # --- 4. Identify Transition and Calculate Energy Difference (ΔE) ---
    # The transition is from the fundamental state (v=0, J=0) to the next
    # accessible state via photon absorption (v=1, J=1).
    # The energy difference is ΔE = E(1,1) - E(0,0) = ħω + ħ²/I.
    vibrational_energy_term = h_bar * omega
    rotational_energy_term = (h_bar ** 2) / I
    delta_E = vibrational_energy_term + rotational_energy_term

    # --- 5. Calculate Photon Momentum (p) ---
    calculated_p = delta_E / c

    # --- 6. Verify the Answer ---
    # The provided final answer is D, which corresponds to p = 1.4 * 10^-28 N*s.
    target_p = 1.4e-28

    # Check if the calculated value is close to the target value.
    # A relative tolerance of 5% is generous enough to account for rounding
    # differences in the options or constants used.
    if math.isclose(calculated_p, target_p, rel_tol=0.05):
        return "Correct"
    else:
        return (
            f"Incorrect. The calculated momentum is p = {calculated_p:.4e} N*s. "
            f"This does not match the value from the selected option D (p = {target_p:.4e} N*s). "
            f"The step-by-step calculation is as follows:\n"
            f"1. Reduced Mass (μ): {mu_kg:.4e} kg\n"
            f"2. Moment of Inertia (I): {I:.4e} kg·m²\n"
            f"3. Transition Energy (ΔE): {delta_E:.4e} J\n"
            f"4. Final Momentum (p = ΔE/c): {calculated_p:.4e} N*s"
        )

# Run the check
result = check_diatomic_molecule_momentum()
print(result)