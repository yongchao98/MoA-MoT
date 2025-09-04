import numpy as np

def check_lienard_wiechert_potentials():
    """
    Checks the correctness of the Liénard-Wiechert potentials from the selected answer.

    It verifies two fundamental constraints:
    1. The potentials must reduce to the standard electrostatic/magnetostatic 
       potentials for a stationary charge (v=0).
    2. The scalar and vector potentials must satisfy the relation A = (v/c^2) * V.

    Returns:
        str: "Correct" if all checks pass, otherwise a string explaining the failure.
    """
    # Define physical constants for the check.
    # The exact values don't matter, but their relationship must be correct.
    q = 1.6e-19  # Charge
    c = 3.0e8      # Speed of light
    epsilon0 = 8.854e-12 # Permittivity of free space
    mu0 = 4 * np.pi * 1e-7 # Permeability of free space

    # Internal sanity check that constants are consistent (c^2 ≈ 1/(μ₀ε₀))
    if not np.isclose(c**2, 1 / (mu0 * epsilon0)):
        return "Error in checker: Physical constants are not consistent."

    # Define the formulas from the proposed answer (Option B) as a function
    def get_potentials_from_answer(charge, velocity, d_vector):
        d_magnitude = np.linalg.norm(d_vector)
        if d_magnitude == 0:
            return float('inf'), np.array([float('inf')] * 3)
        
        d_dot_v = np.dot(d_vector, velocity)
        denominator = d_magnitude * c - d_dot_v
        
        if np.isclose(denominator, 0):
            return float('inf'), np.array([float('inf')] * 3)

        # Scalar Potential V from Option B
        V = (charge * c) / (4 * np.pi * epsilon0 * denominator)
        # Vector Potential A from Option B
        A = (mu0 * charge * c * velocity) / (4 * np.pi * denominator)
        return V, A

    # --- Constraint 1: Check the Static Limit (v = 0) ---
    v_static = np.array([0.0, 0.0, 0.0])
    d_vec_1 = np.array([10.0, 0.0, 0.0]) # An arbitrary distance vector
    d_mag_1 = np.linalg.norm(d_vec_1)

    # Calculate potentials using the answer's formula for a static charge
    V_calc_static, A_calc_static = get_potentials_from_answer(q, v_static, d_vec_1)

    # Calculate the expected potentials for a static charge
    V_expected_static = q / (4 * np.pi * epsilon0 * d_mag_1)
    A_expected_static = np.array([0.0, 0.0, 0.0])

    if not np.isclose(V_calc_static, V_expected_static):
        return (f"Incorrect: The formula for V does not correctly reduce to the Coulomb potential for a stationary charge.\n"
                f"Expected: {V_expected_static}\n"
                f"Got: {V_calc_static}")

    if not np.allclose(A_calc_static, A_expected_static):
        return (f"Incorrect: The formula for A does not correctly reduce to zero for a stationary charge.\n"
                f"Expected: {A_expected_static}\n"
                f"Got: {A_calc_static}")

    # --- Constraint 2: Check the A-V Relationship for a moving charge ---
    v_moving = np.array([0.5 * c, 0.1 * c, -0.2 * c]) # A relativistic velocity
    d_vec_2 = np.array([3.0, -4.0, 5.0]) # An arbitrary distance vector

    # Calculate potentials using the answer's formula for a moving charge
    V_calc_moving, A_calc_moving = get_potentials_from_answer(q, v_moving, d_vec_2)

    # Check if A ≈ (v/c^2) * V
    A_expected_from_V = (v_moving / c**2) * V_calc_moving
    
    if not np.allclose(A_calc_moving, A_expected_from_V, rtol=1e-7):
        return (f"Incorrect: The formulas for V and A do not satisfy the relationship A = (v/c^2) * V.\n"
                f"A calculated directly: {A_calc_moving}\n"
                f"A derived from V:      {A_expected_from_V}")

    # If both constraints are satisfied, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_lienard_wiechert_potentials()
# This will print "Correct" if the answer B is valid, or an error message otherwise.
# print(result)
if result == "Correct":
    # The final output should be just the word "Correct" as per the instructions.
    print("Correct")
else:
    # Or the reason for the failure.
    print(result)