import numpy as np
from scipy.constants import c, epsilon_0, mu_0, pi

def check_correctness_of_answer_B():
    """
    This function programmatically checks the correctness of the formulas for the
    scalar and vector potentials given in Option B.

    The formulas from Option B are:
    V(r,t) = qc / (4πε₀(dc - d·v))
    A(r,t) = μ₀qc v / (4π(dc - d·v))
    """

    # Define a function to compute potentials based on Option B's formulas
    def calculate_potentials_B(q, d_vec, v_vec):
        d = np.linalg.norm(d_vec)
        if d == 0:
            return float('inf'), np.array([float('inf')] * 3)

        d_dot_v = np.dot(d_vec, v_vec)
        denominator = d * c - d_dot_v
        
        if np.isclose(denominator, 0):
            return float('inf'), np.array([float('inf')] * 3)

        V = (q * c) / (4 * pi * epsilon_0 * denominator)
        A_vec = (mu_0 * q * c * v_vec) / (4 * pi * denominator)
        
        return V, A_vec

    # --- Test 1: Static Limit (v = 0) ---
    # The potentials should reduce to the electrostatic case.
    q_test = 1.602e-19  # Elementary charge
    d_test_vec = np.array([1e-9, 2e-9, -2e-9]) # A test distance vector
    v_static = np.array([0.0, 0.0, 0.0])
    
    V_calc_static, A_calc_static = calculate_potentials_B(q_test, d_test_vec, v_static)
    
    # Expected results for a static charge
    d_mag = np.linalg.norm(d_test_vec)
    V_expected_static = q_test / (4 * pi * epsilon_0 * d_mag)
    A_expected_static = np.array([0.0, 0.0, 0.0])
    
    if not np.isclose(V_calc_static, V_expected_static):
        return f"Incorrect. The scalar potential V from Option B fails the static limit test (v=0). Expected {V_expected_static}, but got {V_calc_static}."
    if not np.allclose(A_calc_static, A_expected_static):
        return f"Incorrect. The vector potential A from Option B fails the static limit test (v=0). Expected {A_expected_static}, but got {A_calc_static}."

    # --- Test 2: Internal Consistency (A = v/c^2 * V) ---
    # This relationship must hold for any valid v and d.
    v_general = np.array([0.1*c, -0.2*c, 0.05*c]) # A general relativistic velocity
    
    V_calc_general, A_calc_general = calculate_potentials_B(q_test, d_test_vec, v_general)
    
    # Calculate A independently from V using the consistency relation
    A_from_V_relation = (v_general / c**2) * V_calc_general
    
    # Compare the A calculated from the formula with the A derived from V
    if not np.allclose(A_calc_general, A_from_V_relation, rtol=1e-9):
        return f"Incorrect. Option B fails the internal consistency test. The relation A = (v/c^2)V is not satisfied. A from formula: {A_calc_general}, A derived from V: {A_from_V_relation}."

    # If both fundamental tests pass, the answer is correct.
    # The formulas in B are algebraically equivalent to the standard Liénard-Wiechert potentials.
    return "Correct"

# Run the verification
result = check_correctness_of_answer_B()
print(result)