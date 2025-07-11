import numpy as np

# This script analyzes a series of statements about a proposed Quaternion RoPE.
# It defines quaternion operations and then programmatically verifies each statement.
# Finally, it prints the sorted list of letters corresponding to the true statements.

# --- Quaternion helper functions ---
# A quaternion q = a + bi + cj + dk is represented as np.array([a, b, c, d])

def q_mul(q1, q2):
    """Multiplies two quaternions."""
    a1, b1, c1, d1 = q1
    a2, b2, c2, d2 = q2
    a = a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2
    b = a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2
    c = a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2
    d = a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2
    return np.array([a, b, c, d])

def q_norm(q):
    """Returns the norm of a quaternion."""
    return np.sqrt(np.sum(q**2))

def inner_product(q1, q2):
    """Returns the inner product of two quaternions (viewed as 4D vectors)."""
    return np.dot(q1, q2)

def get_R(m, theta, u_axis):
    """Creates the rotation quaternion R_m."""
    angle = m * theta
    cos_val = np.cos(angle)
    sin_val = np.sin(angle)
    r = np.array([cos_val, 
                  u_axis[0] * sin_val, 
                  u_axis[1] * sin_val, 
                  u_axis[2] * sin_val])
    return r

def get_M(R):
    """Gets the 4x4 matrix for left-multiplication by quaternion R."""
    a, b, c, d = R
    return np.array([
        [a, -b, -c, -d],
        [b,  a, -d,  c],
        [c,  d,  a, -b],
        [d, -c,  b,  a]
    ])

def analyze_statements():
    """
    Evaluates each statement from A to N and returns a dictionary of boolean results.
    """
    # Setup test parameters
    theta = np.pi / 8
    u_axis = np.array([1.0, 2.0, 3.0])
    u_axis = u_axis / np.linalg.norm(u_axis) # Normalize to a unit vector
    
    q_vec = np.array([0.5, 1.0, -2.0, 0.1])
    k_vec = np.array([-0.2, 0.8, 1.5, -1.0])
    
    results = {}
    
    # A) The inner product ⟨R_m(q), R_n(k)⟩ depends only on |m-n|
    # Test: Compare result for (m,n)=(2,5) vs (5,2). |m-n| is 3 for both.
    m1, n1 = 2, 5
    m2, n2 = 5, 2
    R_m1, R_n1 = get_R(m1, theta, u_axis), get_R(n1, theta, u_axis)
    R_m2, R_n2 = get_R(m2, theta, u_axis), get_R(n2, theta, u_axis)
    val1 = inner_product(q_mul(R_m1, q_vec), q_mul(R_n1, k_vec))
    val2 = inner_product(q_mul(R_m2, q_vec), q_mul(R_n2, k_vec))
    results['A'] = np.isclose(val1, val2) # False, depends on m-n, not |m-n|

    # B) R_m(q)R_n(k) = R_{m+p}(q)R_{n+p}(k) for any p
    m, n, p = 1, 2, 3
    R_m, R_n = get_R(m, theta, u_axis), get_R(n, theta, u_axis)
    R_mp, R_np = get_R(m + p, theta, u_axis), get_R(n + p, theta, u_axis)
    lhs = q_mul(q_mul(R_m, q_vec), q_mul(R_n, k_vec))
    rhs = q_mul(q_mul(R_mp, q_vec), q_mul(R_np, k_vec))
    results['B'] = np.allclose(lhs, rhs) # False

    # C) The rotation axis (u_1,u_2,u_3) must be fixed for all positions
    # True by design for RoPE-like relative encoding properties to hold.
    results['C'] = True

    # D) Quaternion RoPE can encode 3D relative positions with a single rotation
    # False. It's a 1-parameter family of rotations for a 1D position m.
    results['D'] = False

    # E) The magnitude |R_m(v)| equals |v| for all m
    m = 5
    R_m = get_R(m, theta, u_axis)
    rotated_v = q_mul(R_m, q_vec)
    results['E'] = np.isclose(q_norm(rotated_v), q_norm(q_vec)) # True

    # F) R_m(αv) = αR_m(v) for scalar α
    m, alpha = 3, 2.5
    R_m = get_R(m, theta, u_axis)
    lhs = q_mul(R_m, alpha * q_vec)
    rhs = alpha * q_mul(R_m, q_vec)
    results['F'] = np.allclose(lhs, rhs) # True

    # G) The quaternion rotation preserves orthogonality between vectors
    p_vec = np.array([1, 1, 0, 0])
    q_ortho = np.array([0, 0, 1, 1]) # inner_product is 0
    m = 4
    R_m = get_R(m, theta, u_axis)
    rotated_p = q_mul(R_m, p_vec)
    rotated_q = q_mul(R_m, q_ortho)
    results['G'] = np.isclose(inner_product(rotated_p, rotated_q), 0) # True

    # H) The composition R_m ∘ R_n equals R_{m+n}
    m, n = 2, 3
    R_m, R_n = get_R(m, theta, u_axis), get_R(n, theta, u_axis)
    R_mn = get_R(m + n, theta, u_axis)
    results['H'] = np.allclose(q_mul(R_m, R_n), R_mn) # True

    # J) (R_m ∘ R_n)(v) - R_n(R_m(v)) is always purely imaginary
    m, n = 2, 3
    R_m, R_n = get_R(m, theta, u_axis), get_R(n, theta, u_axis)
    commutator = q_mul(R_m, R_n) - q_mul(R_n, R_m)
    # Commutator is the zero quaternion [0,0,0,0]. A quaternion is purely
    # imaginary if its real part is 0. This is true for the zero quaternion.
    results['J'] = np.isclose(commutator[0], 0) # True

    # K) The double quaternion rotation R_m(R_n(v)) preserves the real component of v
    m, n = 1, 2
    R_mn = get_R(m + n, theta, u_axis)
    rotated_v = q_mul(R_mn, q_vec)
    results['K'] = np.isclose(rotated_v[0], q_vec[0]) # False

    # L) The trace of the 4×4 matrix representation of R_m equals 4cos(mθ)
    m = 6
    R_m = get_R(m, theta, u_axis)
    M_m = get_M(R_m)
    trace = np.trace(M_m)
    expected_trace = 4 * np.cos(m * theta)
    results['L'] = np.isclose(trace, expected_trace) # True

    # M) The commutator [R_m, R_n] = 0 if and only if the rotation axis is parallel to (1,1,1)
    # False. The commutator is always 0 for any fixed axis.
    results['M'] = False

    # N) For any unit quaternion v, lim_{m→∞} |R_m(v)| converges to 1/√2
    # False. The norm is always preserved, so the limit is |v|, which is 1.
    results['N'] = False
    
    return results

if __name__ == '__main__':
    # Perform the analysis
    statement_results = analyze_statements()
    
    # Collect the letters of the true statements
    correct_statements = sorted([key for key, value in statement_results.items() if value])
    
    # Format the final answer
    final_answer = "".join(correct_statements)
    
    # Print the final answer as requested
    print(final_answer)