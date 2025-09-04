import math

def check_answer():
    """
    This function checks the correctness of the provided answer by recalculating the photon momentum
    required for the specified molecular transition.
    """
    # --- Physical Constants (in SI units) ---
    h_bar = 1.054571817e-34  # Reduced Planck constant in J*s
    c = 2.99792458e8         # Speed of light in m/s
    amu_to_kg = 1.66053906660e-27 # Conversion factor from atomic mass units to kg

    # --- Given Parameters from the Question ---
    Mx_amu = 20.0
    My_amu = 2.0
    R_angstrom = 2.0
    w = 4.0e14  # Angular frequency in rad/s

    # --- Convert given parameters to SI units ---
    Mx_kg = Mx_amu * amu_to_kg
    My_kg = My_amu * amu_to_kg
    R_m = R_angstrom * 1e-10

    # --- Calculation Steps ---

    # The problem asks for the transition from the fundamental state (v=0, J=0)
    # to the next state with the lowest possible energy. As reasoned in the provided text,
    # this corresponds to a rovibrational transition to (v=1, J=1), which is allowed by
    # selection rules (Δv=+1, ΔJ=+1).

    # The energy of this transition is ΔE = E(1,1) - E(0,0).
    # E(v,J) = (v + 1/2)ħω + B*J(J+1)
    # ΔE = [(3/2)ħω + 2B] - [(1/2)ħω] = ħω + 2B

    # 1. Calculate the reduced mass (μ) of the molecule.
    mu = (Mx_kg * My_kg) / (Mx_kg + My_kg)

    # 2. Calculate the moment of inertia (I).
    I = mu * R_m**2

    # 3. Calculate the rotational constant (B).
    B = (h_bar**2) / (2 * I)

    # 4. Calculate the total energy change (ΔE) for the transition.
    delta_E = (h_bar * w) + (2 * B)

    # 5. Calculate the required photon momentum (p = ΔE / c).
    calculated_p = delta_E / c

    # --- Verification ---

    # The provided answer is B) p = 1.4 * 10^(-28) N*s
    answer_p = 1.4e-28

    # Check if the calculated value is close to the answer's value.
    # A relative tolerance of 5% is reasonable for this type of problem.
    if math.isclose(calculated_p, answer_p, rel_tol=0.05):
        return "Correct"
    else:
        # If the answer is incorrect, provide the calculated result and intermediate steps.
        reason = (
            f"The provided answer B (p = {answer_p:.1e} N*s) is incorrect.\n"
            f"The calculation, following the correct physical model of a rovibrational transition from (v=0, J=0) to (v=1, J=1), yields a different result.\n"
            f"Here are the calculation steps:\n"
            f"1. Reduced mass (μ): {mu:.4e} kg\n"
            f"2. Moment of inertia (I): {I:.4e} kg*m^2\n"
            f"3. Rotational constant (B): {B:.4e} J\n"
            f"4. Total energy change (ΔE = ħω + 2B): {delta_E:.4e} J\n"
            f"5. Calculated photon momentum (p = ΔE/c): {calculated_p:.4e} N*s\n"
            f"The calculated momentum {calculated_p:.2e} N*s does not match the value from option B."
        )
        return reason

# Run the check
result = check_answer()
print(result)