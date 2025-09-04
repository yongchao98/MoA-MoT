import math

def check_answer():
    """
    Checks the correctness of the answer for the diatomic molecule momentum problem.
    """
    # --- Define physical constants in SI units ---
    h_bar = 1.054571817e-34  # Reduced Planck constant in J*s
    c = 2.99792458e8         # Speed of light in m/s
    amu_to_kg = 1.660539e-27 # AMU to kg conversion factor

    # --- Given parameters from the question ---
    Mx_amu = 20.0
    My_amu = 2.0
    R_angstrom = 2.0
    omega = 4.0e14 # rad/s

    # --- Convert parameters to SI units ---
    Mx_kg = Mx_amu * amu_to_kg
    My_kg = My_amu * amu_to_kg
    R_m = R_angstrom * 1e-10

    # --- Step 1: Determine the transition ---
    # The problem asks for the transition from the fundamental state to the next state with the lowest possible energy.
    # The fundamental (ground) state is (v=0, J=0).
    # For photon absorption, the selection rules are Δv = +1 and ΔJ = ±1.
    # Starting from J=0, ΔJ=-1 is impossible. So the transition must be ΔJ=+1.
    # Therefore, the transition is from (v=0, J=0) to (v=1, J=1).

    # --- Step 2: Calculate the transition energy (ΔE) ---
    # The energy of a state (v, J) is E(v, J) = ħω(v + 1/2) + (ħ²/2I) * J(J+1).
    # E_initial = E(0, 0) = ħω(1/2)
    # E_final = E(1, 1) = ħω(3/2) + (ħ²/2I) * 1*(1+1) = ħω(3/2) + ħ²/I
    # ΔE = E_final - E_initial = ħω + ħ²/I

    # To calculate ΔE, we first need the reduced mass (μ) and moment of inertia (I).
    
    # --- Step 3: Calculate reduced mass (μ) ---
    mu = (Mx_kg * My_kg) / (Mx_kg + My_kg)

    # --- Step 4: Calculate moment of inertia (I) ---
    I = mu * R_m**2

    # --- Step 5: Calculate the final transition energy (ΔE) ---
    delta_E = (h_bar * omega) + (h_bar**2 / I)

    # --- Step 6: Calculate the photon's momentum (p) ---
    # For a photon, E = pc, so p = E/c.
    p_calculated = delta_E / c

    # --- Step 7: Compare the calculated value with the given answer ---
    # The provided answer is A, which corresponds to p = 1.4 * 10^-28 N*s.
    answer_value = 1.4e-28

    # We check if the calculated value is very close to the answer's value.
    # A tolerance of 5% is reasonable for this kind of problem.
    if abs(p_calculated - answer_value) / answer_value < 0.05:
        return "Correct"
    else:
        return (f"Incorrect. The calculated momentum is approximately {p_calculated:.3e} N*s. "
                f"The value for the chosen option A is {answer_value:.3e} N*s. "
                f"The calculated value does not match the provided answer.")

# Execute the check
result = check_answer()
print(result)