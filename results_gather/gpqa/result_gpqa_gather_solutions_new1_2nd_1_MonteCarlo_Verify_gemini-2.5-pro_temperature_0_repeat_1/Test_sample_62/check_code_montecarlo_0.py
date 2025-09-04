import math

def check_diatomic_momentum():
    """
    This function checks the correctness of the answer to the diatomic molecule problem.
    It calculates the theoretical momentum of a photon required for a specific quantum transition
    and compares it to the value given in the selected answer.
    """

    # --- 1. Define Physical Constants (in SI units) ---
    h_bar = 1.054571817e-34  # Reduced Planck constant in J*s
    c = 2.99792458e8         # Speed of light in m/s
    amu_to_kg = 1.660539e-27 # AMU to kg conversion factor
    angstrom_to_m = 1e-10      # Angstrom to meter conversion factor

    # --- 2. Define Given Parameters from the Question ---
    Mx_amu = 20.0  # Mass of atom X in amu
    My_amu = 2.0   # Mass of atom Y in amu
    R_angstrom = 2.0 # Molecular bond length in angstroms
    omega = 4.0e14   # Angular frequency of vibration in rad/s

    # --- 3. Define the Answer to be Checked ---
    # The final answer provided is 'C', which corresponds to p = 1.4 * 10^-28 N*s
    expected_answer_value = 1.4e-28

    # --- 4. Perform the Physics Calculation ---

    # Convert given parameters to SI units
    Mx_kg = Mx_amu * amu_to_kg
    My_kg = My_amu * amu_to_kg
    R_m = R_angstrom * angstrom_to_m

    # Calculate the reduced mass (mu) of the molecule
    # Formula: mu = (m1 * m2) / (m1 + m2)
    mu = (Mx_kg * My_kg) / (Mx_kg + My_kg)

    # Calculate the moment of inertia (I) of the molecule
    # Formula: I = mu * R^2
    I = mu * R_m**2

    # Calculate the transition energy (Delta_E)
    # The transition is from the ground state (v=0, J=0) to the next allowed state (v=1, J=1)
    # based on selection rules (Δv=+1, ΔJ=+1).
    # E(v,J) = h_bar*omega*(v+1/2) + (h_bar^2 / (2*I)) * J*(J+1)
    # E_initial = E(0,0) = h_bar*omega*(1/2)
    # E_final = E(1,1) = h_bar*omega*(3/2) + (h_bar^2 / I)
    # Delta_E = E_final - E_initial = h_bar*omega + h_bar**2 / I
    delta_E = (h_bar * omega) + (h_bar**2 / I)

    # Calculate the photon's momentum (p)
    # Formula: p = E / c
    p_calculated = delta_E / c

    # --- 5. Verify the Correctness of the Answer ---
    # We check if the calculated momentum is very close to the expected answer's value.
    # A small tolerance (e.g., 2%) is used to account for potential rounding in the problem's options.
    tolerance = 0.02
    if abs(p_calculated - expected_answer_value) / expected_answer_value < tolerance:
        return "Correct"
    else:
        error_message = (
            f"Incorrect. The provided answer corresponds to a momentum of {expected_answer_value:.3e} N*s. "
            f"However, the calculation based on the problem's parameters yields a momentum of {p_calculated:.3e} N*s.\n"
            f"The calculation steps are as follows:\n"
            f"Reduced Mass (μ) = {mu:.4e} kg\n"
            f"Moment of Inertia (I) = {I:.4e} kg*m^2\n"
            f"Transition Energy (ΔE) = {delta_E:.4e} J\n"
            f"Calculated Momentum (p) = {p_calculated:.4e} N*s"
        )
        return error_message

# Run the check and print the result
result = check_diatomic_momentum()
print(result)