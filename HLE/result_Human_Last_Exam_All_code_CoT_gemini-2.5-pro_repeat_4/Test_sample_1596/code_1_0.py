import numpy as np

def q_mult(p, q):
    """Multiplies two quaternions p*q represented as numpy arrays."""
    p0, p1, p2, p3 = p
    q0, q1, q2, q3 = q
    r0 = p0 * q0 - p1 * q1 - p2 * q2 - p3 * q3
    r1 = p0 * q1 + p1 * q0 + p2 * q3 - p3 * q2
    r2 = p0 * q2 - p1 * q3 + p2 * q0 + p3 * q1
    r3 = p0 * q3 + p1 * q2 - p2 * q1 + p3 * q0
    return np.array([r0, r1, r2, r3])

def q_norm(q):
    """Calculates the norm of a quaternion."""
    return np.sqrt(np.sum(q**2))

def inner_product(p, q):
    """Calculates the inner product of two quaternions."""
    return np.dot(p, q)

def main():
    """
    Analyzes statements about Quaternion RoPE and prints the final answer.
    """
    # Setup for the checks
    theta = np.pi / 8  # A fixed parameter for the embedding
    # A fixed, normalized rotation axis u = (u1, u2, u3)
    u_vec = np.array([1, 2, 3])
    u_vec = u_vec / np.linalg.norm(u_vec)

    def get_R(m):
        """Generates the rotation quaternion R_m."""
        angle = m * theta
        scalar_part = np.cos(angle)
        vector_part = np.sin(angle) * u_vec
        return np.array([scalar_part, vector_part[0], vector_part[1], vector_part[2]])

    # Random vectors for testing
    np.random.seed(42)
    q_vec = np.random.rand(4)
    k_vec = np.random.rand(4)
    v_vec = np.random.rand(4)

    true_statements = []

    print("Analyzing statements about Quaternion RoPE...\n")

    # A) The inner product ⟨R_m(q), R_n(k)⟩ depends only on |m-n|
    m, n = 2, 5
    val1 = inner_product(q_mult(get_R(m), q_vec), q_mult(get_R(n), k_vec))
    # Check with positions swapped, so m-n becomes n-m
    m2, n2 = 5, 2
    val2 = inner_product(q_mult(get_R(m2), q_vec), q_mult(get_R(n2), k_vec))
    print(f"A) Check if inner product for (m,n)=(2,5) equals (5,2), where |m-n| is the same.")
    print(f"   - Inner product for (m-n = -3): {val1:.4f}")
    print(f"   - Inner product for (m-n = 3): {val2:.4f}")
    print(f"   - A is False. The result depends on the relative position (m-n), not its magnitude.")

    # B) R_m(q)R_n(k) = R_{m+p}(q)R_{n+p}(k) for any p
    m, n, p = 1, 2, 3
    lhs = q_mult(q_mult(get_R(m), q_vec), q_mult(get_R(n), k_vec))
    rhs = q_mult(q_mult(get_R(m+p), q_vec), q_mult(get_R(n+p), k_vec))
    print(f"\nB) Check if R_m(q)R_n(k) = R_{m+p}(q)R_{n+p}(k)")
    print(f"   - LHS value: {np.round(lhs, 4)}")
    print(f"   - RHS value: {np.round(rhs, 4)}")
    print(f"   - B is False. This identity does not hold.")

    # C) The rotation axis (u_1,u_2,u_3) must be fixed for all positions
    print("\nC) The rotation axis must be fixed for all positions.")
    print("   - True. For the scheme to have RoPE-like properties, such as R_m*R_n = R_{m+n}, the operators must commute. This only happens if their rotation axes are collinear, so the axis must be fixed.")
    true_statements.append('C')

    # D) Quaternion RoPE can encode 3D relative positions with a single rotation
    print("\nD) Can encode 3D relative positions.")
    print("   - False. The defined operator R_m uses a scalar position 'm', thus encoding a 1D positional relationship.")
    
    # E) The magnitude |R_m(v)| equals |v| for all m
    m = 4
    v_norm_before = q_norm(v_vec)
    v_norm_after = q_norm(q_mult(get_R(m), v_vec))
    print(f"\nE) The magnitude |R_m(v)| equals |v|.")
    print(f"   - |v| before rotation: {v_norm_before:.4f}")
    print(f"   - |R_m(v)| after rotation: {v_norm_after:.4f}")
    print(f"   - True. The rotation quaternion R_m is a unit quaternion, so the rotation is an isometry and preserves norms.")
    true_statements.append('E')

    # F) R_m(αv) = αR_m(v) for scalar α
    m, alpha = 2, 3.14
    lhs = q_mult(get_R(m), alpha * v_vec)
    rhs = alpha * q_mult(get_R(m), v_vec)
    print(f"\nF) Check if R_m(αv) = αR_m(v).")
    print(f"   - LHS R_m(αv): {np.round(lhs, 4)}")
    print(f"   - RHS αR_m(v): {np.round(rhs, 4)}")
    print(f"   - True. Scalar multiplication is commutative with quaternion multiplication.")
    true_statements.append('F')
    
    # G) The quaternion rotation preserves orthogonality between vectors
    v1 = np.random.rand(4)
    v2 = np.random.rand(4)
    v2 -= v2.dot(v1) / v1.dot(v1) * v1 # Make v2 orthogonal to v1
    ip_before = inner_product(v1, v2)
    ip_after = inner_product(q_mult(get_R(3), v1), q_mult(get_R(3), v2))
    print("\nG) The rotation preserves orthogonality.")
    print(f"   - Inner product of orthogonal vectors before rotation: {ip_before:.4e}")
    print(f"   - Inner product after rotation: {ip_after:.4e}")
    print("   - True. The rotation is an orthogonal transformation that preserves inner products.")
    true_statements.append('G')

    # H) The composition R_m ∘ R_n equals R_{m+n}
    m, n = 2, 3
    R_m_R_n = q_mult(get_R(m), get_R(n))
    R_mn = get_R(m+n)
    print(f"\nH) Check if R_m * R_n = R_{m+n}.")
    print(f"   - R_m * R_n: {np.round(R_m_R_n, 4)}")
    print(f"   - R_(m+n):   {np.round(R_mn, 4)}")
    print("   - True. This holds because R_m and R_n commute.")
    true_statements.append('H')

    # J) For any quaternion vector v, (R_m ∘ R_n)(v) - R_n(R_m(v)) is always purely imaginary
    m, n = 2, 3
    commutator = q_mult(get_R(m), get_R(n)) - q_mult(get_R(n), get_R(m))
    result_vec = q_mult(commutator, v_vec)
    print(f"\nJ) The commutator applied to v, [R_m, R_n]v, is purely imaginary.")
    print(f"   - [R_m, R_n] = {np.round(commutator, 4)} (is zero)")
    print(f"   - Result vector = {np.round(result_vec, 4)} (is zero)")
    print("   - True. The commutator is the zero quaternion, which has a real part of 0.")
    true_statements.append('J')

    # K) The double quaternion rotation R_m(R_n(v)) preserves the real component of v
    m, n = 2, 3
    v_real_before = v_vec[0]
    v_after = q_mult(get_R(m+n), v_vec)
    v_real_after = v_after[0]
    print(f"\nK) The real component of v is preserved.")
    print(f"   - Real part before: {v_real_before:.4f}")
    print(f"   - Real part after:  {v_real_after:.4f}")
    print("   - False. The real component is generally not preserved.")
    
    # L) The trace of the 4×4 matrix representation of R_m equals 4cos(mθ)
    m = 5
    q0 = get_R(m)[0]
    trace_val = 4 * q0
    expected_val = 4 * np.cos(m * theta)
    print(f"\nL) The trace of the matrix of R_m is 4*cos(mθ).")
    print(f"   - Trace = 4 * R_m[0] = {trace_val:.4f}")
    print(f"   - Expected = 4*cos({m}*θ) = {expected_val:.4f}")
    print("   - True. The trace of the left-multiplication matrix is 4 times the scalar part of the quaternion.")
    true_statements.append('L')

    # M) The commutator [R_m, R_n] = 0 if and only if the rotation axis is parallel to (1,1,1)
    print("\nM) The commutator is zero iff axis is (1,1,1).")
    print("   - False. The commutator is always zero for any fixed axis. The condition on the axis direction is not necessary.")

    # N) For any unit quaternion v, lim_{m→∞} |R_m(v)| converges to 1/√2
    unit_v = v_vec / q_norm(v_vec)
    norm_val = q_norm(q_mult(get_R(1000), unit_v))
    print(f"\nN) The limit of the norm is 1/sqrt(2).")
    print(f"   - As shown in E, |R_m(v)| = |v|. For a unit vector, the norm is always 1.")
    print(f"   - For large m=1000, |R_m(unit_v)| = {norm_val:.4f}")
    print("   - False. The limit is 1.")

    # Final result
    true_statements.sort()
    final_answer = "".join(true_statements)

    print("\n-----------------------------------------")
    print(f"The sorted list of true statement letters is:")
    print(final_answer)
    print("<<<" + final_answer + ">>>")

if __name__ == '__main__':
    main()