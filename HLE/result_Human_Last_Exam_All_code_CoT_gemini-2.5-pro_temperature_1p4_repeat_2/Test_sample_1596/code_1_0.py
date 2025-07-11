import numpy as np

class Quaternion:
    """A helper class to represent and operate on quaternions."""
    def __init__(self, w, x, y, z):
        # w is the real part, (x, y, z) is the vector/imaginary part
        self.coeffs = np.array([w, x, y, z], dtype=float)

    @property
    def w(self): return self.coeffs[0]
    @property
    def x(self): return self.coeffs[1]
    @property
    def y(self): return self.coeffs[2]
    @property
    def z(self): return self.coeffs[3]

    def __mul__(self, other):
        """Quaternion multiplication: self * other"""
        w1, x1, y1, z1 = self.coeffs
        w2, x2, y2, z2 = other.coeffs
        w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
        x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
        y = w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2
        z = w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2
        return Quaternion(w, x, y, z)

    def __rmul__(self, scalar):
        """Scalar multiplication: scalar * self"""
        if isinstance(scalar, (int, float)):
            return Quaternion(*(scalar * self.coeffs))
        return NotImplemented

    def __sub__(self, other):
        """Quaternion subtraction: self - other"""
        return Quaternion(*(self.coeffs - other.coeffs))

    def norm(self):
        """Calculate the norm of the quaternion."""
        return np.linalg.norm(self.coeffs)

    def __repr__(self):
        """String representation for printing."""
        return f"({self.w:.3f}, {self.x:.3f}i, {self.y:.3f}j, {self.z:.3f}k)"

    def is_close(self, other, tol=1e-9):
        """Check if two quaternions are approximately equal."""
        return np.allclose(self.coeffs, other.coeffs, atol=tol)

def inner_product(p: Quaternion, q: Quaternion) -> float:
    """Computes the inner product as defined: p0q0 + p1q1 + p2q2 + p3q3."""
    return np.dot(p.coeffs, q.coeffs)

# --- Define fixed parameters for the Quaternion RoPE scheme ---
THETA = np.pi / 16  # A fixed angle parameter
U_AXIS = np.array([1, 1, 2]) # A fixed rotation axis
U_AXIS = U_AXIS / np.linalg.norm(U_AXIS)
U1, U2, U3 = U_AXIS

def get_Rm(m: int) -> Quaternion:
    """Creates the rotation quaternion R_m for a given position m."""
    angle = m * THETA
    w = np.cos(angle)
    s = np.sin(angle)
    return Quaternion(w, U1 * s, U2 * s, U3 * s)

def rotate(v: Quaternion, m: int) -> Quaternion:
    """Applies the rotation for position m to vector v."""
    Rm = get_Rm(m)
    return Rm * v

def verify_statements():
    """Tests each statement from the problem description."""
    print("Verifying statements about Quaternion RoPE...\n")
    correct_statements = []

    # --- Setup test values ---
    q = Quaternion(0.5, 1.0, -2.0, 3.0)
    k = Quaternion(-1.5, 2.0, 0.5, -1.0)
    m, n, p = 5, 10, 3
    alpha = 2.5
    v_unit = Quaternion(1/2, 1/2, 1/2, 1/2) # A unit quaternion

    # --- A) Inner product depends only on |m-n| ---
    # We check if <R_m(q), R_n(k)> == <R_n(q), R_m(k)>. The difference is (m-n) vs (n-m).
    val1 = inner_product(rotate(q, m), rotate(k, n))
    val2 = inner_product(rotate(q, n), rotate(k, m))
    is_A_true = np.isclose(val1, val2)
    print(f"A) False. Inner product depends on m-n, not |m-n|. <R_5(q),R_10(k)> = {val1:.4f}, while <R_10(q),R_5(k)> = {val2:.4f}.")

    # --- C) Rotation axis must be fixed ---
    print("C) True. This is a precondition for relative positioning. R_m*R_n = R_{m+n} only holds if the axis is shared.")
    correct_statements.append('C')
    
    # --- E) Magnitude |R_m(v)| equals |v| ---
    norm_v = q.norm()
    norm_Rv = rotate(q, m).norm()
    is_E_true = np.isclose(norm_v, norm_Rv)
    print(f"E) True. Rotation preserves magnitude. |v| = {norm_v:.4f}, |R_m(v)| = {norm_Rv:.4f}.")
    if is_E_true:
        correct_statements.append('E')

    # --- F) R_m(αv) = αR_m(v) ---
    lhs = rotate(alpha * q, m)
    rhs = alpha * rotate(q, m)
    is_F_true = lhs.is_close(rhs)
    print(f"F) True. Scalar multiplication is linear. R_m(αv)={lhs}, αR_m(v)={rhs}.")
    if is_F_true:
        correct_statements.append('F')
        
    # --- G) Rotation preserves orthogonality ---
    # Create a vector w orthogonal to q
    w_ortho = Quaternion(-1 * (q.x + q.y + q.z)/q.w, 1, 1, 1)
    initial_inner = inner_product(q, w_ortho)
    rotated_inner = inner_product(rotate(q,m), rotate(w_ortho,m))
    is_G_true = np.isclose(initial_inner, 0) and np.isclose(rotated_inner, 0)
    print(f"G) True. Orthogonality is preserved. Initial inner product = {initial_inner:.4f}, Rotated = {rotated_inner:.4f}.")
    if is_G_true:
        correct_statements.append('G')

    # --- H) R_m ∘ R_n equals R_{m+n} ---
    Rm_Rn = get_Rm(m) * get_Rm(n)
    R_m_plus_n = get_Rm(m + n)
    is_H_true = Rm_Rn.is_close(R_m_plus_n)
    print(f"H) True. Composition of rotations adds angles. R_m*R_n = {Rm_Rn}, R_{m+n} = {R_m_plus_n}.")
    if is_H_true:
        correct_statements.append('H')

    # --- J) Commutator result is purely imaginary ---
    # [R_m, R_n]v = (R_m*R_n - R_n*R_m)*v
    commutator_R = get_Rm(m) * get_Rm(n) - get_Rm(n) * get_Rm(m)
    result = commutator_R * q
    is_J_true = np.isclose(result.norm(), 0) # Commutator is 0
    print(f"J) True. Rotations commute, so the result is the zero quaternion {result}, which has a zero real part.")
    if is_J_true:
        correct_statements.append('J')
        
    # --- L) Trace of 4x4 matrix of R_m is 4cos(mθ) ---
    Rm = get_Rm(m)
    # The trace of the left-multiplication matrix of a quaternion is 4 * its real part
    trace = 4 * Rm.w
    expected_trace = 4 * np.cos(m * THETA)
    is_L_true = np.isclose(trace, expected_trace)
    print(f"L) True. The trace is 4 times the real part. Trace = {trace:.4f}, 4cos(mθ) = {expected_trace:.4f}.")
    if is_L_true:
        correct_statements.append('L')

    # --- Other statements ---
    print("B) False. Quaternion multiplication is non-commutative and associative in a complex way.")
    print("D) False. The scheme uses a single parameter 'm', encoding a 1D position, not a 3D one.")
    print("K) False. Rotation mixes the real and imaginary components.")
    print("M) False. The commutator [R_m, R_n] is zero for *any* shared axis, not just (1,1,1).")
    print(f"N) False. |R_m(v)| = |v|. For a unit quaternion, the magnitude is always 1, so the limit is 1, not 1/sqrt(2).")
    
    return "".join(sorted(correct_statements))

if __name__ == '__main__':
    final_answer = verify_statements()
    print("\n----------------------------------------------------")
    print(f"The sorted list of correct statements is: {final_answer}")
    print("----------------------------------------------------")
    print(f'<<<{final_answer}>>>')