import numpy as np

def solve_relativity_problem():
    """
    Uses a computational approach to find the relative speed and total energy
    for the given relativistic system, then verifies against the options.
    This fulfills the prompt's requirement for sampling/exploration followed by
    deterministic verification.
    """
    # --- (a) Sample / Explore: Calculate the physical quantities ---
    # This is a deterministic calculation, which we'll treat as our "exploration"
    # of the solution space.

    # Given parameters (in units of m, c)
    m1_rest = 2.0
    m2_rest = 3.0
    v1 = 0.6  # Astronaut 1's speed in units of c
    v2 = 0.5  # Astronaut 2's speed in units of c
    c = 1.0   # For simplicity in calculation

    # 1. Calculate relative speed (v_rel)
    # Using the relativistic velocity subtraction formula. The absolute value gives the speed.
    v_rel_calc = abs((v1 - v2) / (1 - (v1 * v2) / c**2))

    # 2. Calculate total energy (E_total)
    # Energy E = gamma * m_rest * c^2, where gamma = 1 / sqrt(1 - v^2/c^2)

    # Astronaut 1
    gamma1 = 1 / np.sqrt(1 - v1**2 / c**2)
    E1 = gamma1 * m1_rest * c**2

    # Astronaut 2
    gamma2 = 1 / np.sqrt(1 - v2**2 / c**2)
    E2 = gamma2 * m2_rest * c**2

    # System's Total Energy (in units of mc^2)
    E_total_calc = E1 + E2

    # --- (b) Narrow candidates: Define the options ---
    candidates = {
        "A": {"v_rel": 0.14, "E": 5.96},
        "B": {"v_rel": 1.1, "E": 1.0},
        "C": {"v_rel": 0.1, "E": 4.96},
        "D": {"v_rel": 0.14, "E": 5.0},
    }

    # --- (c) Run exact checks ---
    verified_answer = None
    # Use a tolerance for floating point comparison due to rounding in options
    tolerance = 1e-2

    for option, values in candidates.items():
        # Check if calculated values are close to the option's values
        v_rel_match = np.isclose(v_rel_calc, values["v_rel"], atol=tolerance)
        E_total_match = np.isclose(E_total_calc, values["E"], atol=tolerance)

        if v_rel_match and E_total_match:
            verified_answer = option
            break

    # --- (d) Print the final verified answer ---
    print("--- Calculation Results ---")
    print(f"Calculated relative speed (v_rel/c): {v_rel_calc:.4f}")
    print(f"Calculated total energy (E / mc^2): {E_total_calc:.4f}")
    print("\n--- Verification Against Options ---")

    if verified_answer:
        print(f"Option {verified_answer} matches the calculations within a tolerance of {tolerance}.")
        print(f"Final Verified Answer: v_rel ≈ {candidates[verified_answer]['v_rel']:.2f}c, E ≈ {candidates[verified_answer]['E']:.2f} mc^2")
    else:
        print("No candidate option passed the verification check.")

# Run the solver
solve_relativity_problem()