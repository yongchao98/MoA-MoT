import math

def check_answer():
    """
    This function checks the correctness of the provided answer to the physics problem.
    It recalculates the photon momentum based on the given parameters and compares it
    to the value in the selected option.
    """

    # --- Physical Constants (in SI units) ---
    # Reduced Planck constant (J*s)
    h_bar = 1.054571817e-34
    # Speed of light (m/s)
    c = 2.99792458e8
    # Atomic mass unit to kg conversion factor
    amu_to_kg = 1.66053906660e-27

    # --- Given Parameters from the Question ---
    # Mass of atom X (amu)
    Mx_amu = 20.0
    # Mass of atom Y (amu)
    My_amu = 2.0
    # Molecular bond length (angstroms)
    R_angstrom = 2.0
    # Angular frequency of vibration (rad/s)
    omega = 4.0e14

    # --- Convert inputs to SI units ---
    Mx_kg = Mx_amu * amu_to_kg
    My_kg = My_amu * amu_to_kg
    R_m = R_angstrom * 1e-10

    # --- Theoretical Model Verification ---
    # The logic for the transition is as follows:
    # Initial state (fundamental): v=0, J=0.
    # Selection rules for photon absorption: Δv=+1, ΔJ=±1.
    # Since J_initial=0, only ΔJ=+1 is possible.
    # Final state (next lowest energy): v=1, J=1.
    # Energy of transition ΔE = E(1,1) - E(0,0).
    # E(v,J) = (v + 1/2)ħω + B*J(J+1), where B = ħ²/(2I).
    # ΔE = [(3/2)ħω + 2B] - [(1/2)ħω] = ħω + 2B.
    # Since 2B = ħ²/I, the formula for ΔE is: ΔE = ħω + ħ²/I.
    # The theoretical steps in the provided answer are correct.

    # --- Step-by-step Calculation ---

    # 1. Calculate the reduced mass (μ) in kg
    # μ = (Mx * My) / (Mx + My)
    mu = (Mx_kg * My_kg) / (Mx_kg + My_kg)

    # 2. Calculate the moment of inertia (I) in kg*m^2
    # I = μ * R^2
    I = mu * R_m**2

    # 3. Calculate the energy difference (ΔE) in Joules
    vibrational_energy_term = h_bar * omega
    rotational_energy_term = h_bar**2 / I
    delta_E = vibrational_energy_term + rotational_energy_term

    # 4. Calculate the photon's momentum (p) in N*s
    # p = ΔE / c
    p_calculated = delta_E / c

    # --- Verification ---
    # The provided answer is B) p = 1.4 * 10^(-28) N*s
    answer_value = 1.4e-28

    # Check if the calculated momentum is close to the answer's value.
    # A relative tolerance of 5% is reasonable for multiple-choice physics problems.
    if math.isclose(p_calculated, answer_value, rel_tol=0.05):
        return "Correct"
    else:
        return (f"Incorrect. The calculated momentum is p = {p_calculated:.3e} N*s, "
                f"which is not close to the value in the selected answer B) {answer_value:.1e} N*s. "
                f"The provided solution's reasoning is correct, but the final numerical result does not match.")

# The code will return "Correct" if the calculation matches the answer.
# Otherwise, it will return a reason for the discrepancy.
result = check_answer()
# print(result)
# Expected output: Correct