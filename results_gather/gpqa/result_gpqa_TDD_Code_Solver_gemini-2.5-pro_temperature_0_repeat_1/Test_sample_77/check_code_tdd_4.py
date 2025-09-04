import numpy as np
from scipy.constants import c, epsilon_0, mu_0, pi

def check_answer():
    """
    Checks the correctness of the formulas in option B for the Liénard-Wiechert potentials.
    The function runs three tests:
    1. Stationary Charge: For v=0, potentials must reduce to the static case (Coulomb potential and zero vector potential).
    2. Perpendicular Motion: For d perpendicular to v, the formulas simplify, and A must be non-zero.
    3. Internal Consistency: The potentials must satisfy the relation A = (v/c^2) * V.
    """

    # Formulas from option B
    def calculate_potentials_B(q, d_vec, v_vec):
        """
        Implements the formulas from option B.
        V(r,t) = qc / (4*pi*epsilon_0 * (d*c - d.v))
        A(r,t) = mu_0 * q * c * v / (4*pi * (d*c - d.v))
        """
        d = np.linalg.norm(d_vec)
        if d == 0:
            return float('inf'), np.array([float('inf')] * len(v_vec))

        d_dot_v = np.dot(d_vec, v_vec)
        
        # Denominator for V
        denominator_V = 4 * pi * epsilon_0 * (d * c - d_dot_v)
        if np.isclose(denominator_V, 0):
            V = float('inf')
        else:
            V = (q * c) / denominator_V

        # Denominator for A
        denominator_A = 4 * pi * (d * c - d_dot_v)
        if np.isclose(denominator_A, 0):
            A_vec = np.array([float('inf')] * len(v_vec))
        else:
            A_vec = (mu_0 * q * c * v_vec) / denominator_A
            
        return V, A_vec

    # --- Test Setup ---
    q_test = 1.602e-19  # Elementary charge in Coulombs

    # --- Test 1: Stationary Charge (v=0) ---
    d_vec_stat = np.array([1e-9, 0.0, 0.0])
    v_vec_stat = np.array([0.0, 0.0, 0.0])
    V_calc_stat, A_calc_stat = calculate_potentials_B(q_test, d_vec_stat, v_vec_stat)
    
    # Expected values for stationary charge
    d_stat = np.linalg.norm(d_vec_stat)
    expected_V_stat = q_test / (4 * pi * epsilon_0 * d_stat)
    expected_A_stat = np.zeros(3)
    
    if not np.isclose(V_calc_stat, expected_V_stat):
        return f"Incorrect. Test 'Stationary Charge' failed for V. For v=0, expected V={expected_V_stat}, but got {V_calc_stat}."
    if not np.allclose(A_calc_stat, expected_A_stat):
        return f"Incorrect. Test 'Stationary Charge' failed for A. For v=0, expected A={expected_A_stat}, but got {A_calc_stat}."

    # --- Test 2: Perpendicular Motion (d . v = 0) ---
    d_vec_perp = np.array([2e-9, 0.0, 0.0])
    v_vec_perp = np.array([0.0, 0.5 * c, 0.0])
    V_calc_perp, A_calc_perp = calculate_potentials_B(q_test, d_vec_perp, v_vec_perp)
    
    # Expected values for perpendicular motion
    d_perp = np.linalg.norm(d_vec_perp)
    # Denominator simplifies: d*c - d.v = d*c.
    # V = qc / (4*pi*eps0*dc) = q / (4*pi*eps0*d)
    expected_V_perp = q_test / (4 * pi * epsilon_0 * d_perp)
    # A = mu0*qc*v / (4*pi*dc) = mu0*q*v / (4*pi*d)
    expected_A_perp = (mu_0 * q_test * v_vec_perp) / (4 * pi * d_perp)
    
    if not np.isclose(V_calc_perp, expected_V_perp):
        return f"Incorrect. Test 'Perpendicular Motion' failed for V. For d.v=0, expected V={expected_V_perp}, but got {V_calc_perp}."
    if np.allclose(A_calc_perp, np.zeros(3)):
        return f"Incorrect. Test 'Perpendicular Motion' failed for A. A should be non-zero but was calculated as zero."
    if not np.allclose(A_calc_perp, expected_A_perp):
        return f"Incorrect. Test 'Perpendicular Motion' failed for A. Expected A={expected_A_perp}, but got {A_calc_perp}."

    # --- Test 3: Internal Consistency (A = v/c^2 * V) ---
    d_vec_gen = np.array([3e-9, -4e-9, 1e-9])
    v_vec_gen = np.array([0.1*c, 0.2*c, -0.05*c])
    V_calc_gen, A_calc_gen = calculate_potentials_B(q_test, d_vec_gen, v_vec_gen)
    
    # Expected A based on the consistency relation
    expected_A_gen = (v_vec_gen / c**2) * V_calc_gen
    
    if not np.allclose(A_calc_gen, expected_A_gen, rtol=1e-9):
        return f"Incorrect. Test 'Internal Consistency' failed. The relation A = (v/c^2)*V is not satisfied. Calculated A={A_calc_gen}, but expected A based on V is {expected_A_gen}."

    return "Correct"

# Run the check
result = check_answer()
print(result)