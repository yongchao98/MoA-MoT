import math
import scipy.constants as const

def check_correctness():
    """
    Checks the correctness of the answer to the diatomic molecule problem.
    """
    # --- Given Parameters ---
    Mx_amu = 20.0
    My_amu = 2.0
    R_m = 2.0e-10  # 2 angstroms in meters
    omega_rad_s = 4.0e14

    # --- Physical Constants from SciPy ---
    h_bar_Js = const.hbar
    c_ms = const.c
    amu_to_kg = const.physical_constants['atomic mass unit-kilogram relationship'][0]

    # --- Step 1: Calculate Reduced Mass (μ) in kg ---
    # The formula for reduced mass is μ = (m1 * m2) / (m1 + m2)
    mu_amu = (Mx_amu * My_amu) / (Mx_amu + My_amu)
    mu_kg = mu_amu * amu_to_kg

    # --- Step 2: Calculate Moment of Inertia (I) in kg*m^2 ---
    # The formula for moment of inertia is I = μ * R^2
    I_kg_m2 = mu_kg * R_m**2

    # --- Step 3: Calculate the Transition Energy (ΔE) in Joules ---
    # The problem describes a transition from the fundamental state (v=0, J=0)
    # to the next lowest accessible state via photon absorption.
    # The selection rules are Δv=+1 and ΔJ=±1.
    # From J=0, only ΔJ=+1 is possible. So the final state is (v=1, J=1).
    # The energy of a state is E(v,J) = ħω(v+1/2) + (ħ²/2I)J(J+1).
    # E_initial = E(0,0) = (1/2)ħω
    # E_final = E(1,1) = (3/2)ħω + (ħ²/2I)*1*(2) = (3/2)ħω + ħ²/I
    # ΔE = E_final - E_initial = ħω + ħ²/I
    vibrational_energy_term = h_bar_Js * omega_rad_s
    rotational_energy_term = h_bar_Js**2 / I_kg_m2
    delta_E_J = vibrational_energy_term + rotational_energy_term

    # --- Step 4: Calculate the Photon's Momentum (p) in N*s ---
    # The energy-momentum relation for a photon is E = pc, so p = E/c.
    calculated_p_Ns = delta_E_J / c_ms

    # --- Step 5: Compare with the provided answer ---
    # The provided answer is D, which corresponds to p = 1.4 * 10^-28 N*s.
    expected_p_Ns = 1.4e-28

    # We check if the calculated value is close to the expected value.
    # A relative tolerance of 2% is reasonable to account for rounding in the option.
    if math.isclose(calculated_p_Ns, expected_p_Ns, rel_tol=0.02):
        return "Correct"
    else:
        return (f"Incorrect. The calculated momentum is {calculated_p_Ns:.4e} N*s, "
                f"which does not match the value from option D ({expected_p_Ns:.1e} N*s).")

# Run the check
result = check_correctness()
print(result)