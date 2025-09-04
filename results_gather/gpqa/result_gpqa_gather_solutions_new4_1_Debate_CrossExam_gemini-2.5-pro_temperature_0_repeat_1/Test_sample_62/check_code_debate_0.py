import math

def check_correctness():
    """
    This function checks the correctness of the selected answer by recalculating the photon momentum
    based on the provided physical parameters and the quantum rotation-vibration model.
    """
    # --- Physical Constants (using high precision values) ---
    h_bar = 1.054571817e-34      # Reduced Planck constant (J·s)
    amu_to_kg = 1.66053906660e-27 # Atomic mass unit to kg conversion
    c = 299792458.0              # Speed of light (m/s)
    angstrom_to_m = 1e-10        # Angstrom to meter conversion

    # --- Given Parameters from the Question ---
    Mx_amu = 20.0
    My_amu = 2.0
    R_angstrom = 2.0
    omega = 4.0e14  # rad/s

    # --- The answer to check is B, which corresponds to p = 1.4e-28 N*s ---
    expected_p = 1.4e-28

    # --- Step 1: Calculate Reduced Mass (μ) in SI units ---
    Mx_kg = Mx_amu * amu_to_kg
    My_kg = My_amu * amu_to_kg
    mu = (Mx_kg * My_kg) / (Mx_kg + My_kg)

    # --- Step 2: Calculate Moment of Inertia (I) in SI units ---
    R_m = R_angstrom * angstrom_to_m
    I = mu * (R_m ** 2)

    # --- Step 3: Calculate the Transition Energy (ΔE) ---
    # The transition is from the ground state (v=0, J=0) to the first excited
    # rovibrational state (v=1, J=1).
    # The energy difference is ΔE = E(1,1) - E(0,0) = ħω + ħ²/I.
    delta_E_vibrational = h_bar * omega
    delta_E_rotational = (h_bar ** 2) / I
    delta_E = delta_E_vibrational + delta_E_rotational

    # --- Step 4: Calculate the Photon's Momentum (p) ---
    calculated_p = delta_E / c

    # --- Step 5: Compare the calculated value with the expected answer ---
    # We use a relative tolerance to account for floating-point inaccuracies
    # and potential rounding in the options. A 5% tolerance is reasonable.
    if math.isclose(calculated_p, expected_p, rel_tol=0.05):
        return "Correct"
    else:
        return (f"Incorrect. The calculation based on the rovibrational transition "
                f"(v=0, J=0) -> (v=1, J=1) yields a momentum of p = {calculated_p:.3e} N·s. "
                f"The selected answer B corresponds to p = {expected_p:.3e} N·s. "
                f"The calculated value does not match the selected option within a 5% tolerance.")

# Execute the check
result = check_correctness()
print(result)