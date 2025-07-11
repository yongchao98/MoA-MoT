import numpy as np

# Step 1: Implement a helper class for Quaternion mathematics
class Quaternion:
    """A simple Quaternion class for the verification task."""
    def __init__(self, w, x, y, z):
        self.w = float(w) # real part
        self.x = float(x) # i part
        self.y = float(y) # j part
        self.z = float(z) # k part
        self.coeffs = np.array([self.w, self.x, self.y, self.z])

    def __mul__(self, other):
        """Quaternion multiplication."""
        w1, x1, y1, z1 = self.coeffs
        w2, x2, y2, z2 = other.coeffs
        w = w1*w2 - x1*x2 - y1*y2 - z1*z2
        x = w1*x2 + x1*w2 + y1*z2 - z1*y2
        y = w1*y2 - x1*z2 + y1*w2 + z1*x2
        z = w1*z2 + x1*y2 - y1*x2 + z1*w2
        return Quaternion(w, x, y, z)
    
    def __rmul__(self, scalar):
        """Scalar multiplication."""
        if isinstance(scalar, (int, float)):
            return Quaternion(scalar * self.w, scalar * self.x, scalar * self.y, scalar * self.z)
        raise TypeError(f"Unsupported operand type for scalar multiplication.")

    def __sub__(self, other):
        """Quaternion subtraction."""
        return Quaternion(self.w - other.w, self.x - other.x, self.y - other.y, self.z - other.z)

    @property
    def norm(self):
        """Quaternion norm."""
        return np.linalg.norm(self.coeffs)

    def is_close(self, other, tol=1e-9):
        """Check if two quaternions are close."""
        return (self - other).norm < tol

def inner_product(q1, q2):
    """Euclidean inner product for quaternions as 4D vectors."""
    return np.dot(q1.coeffs, q2.coeffs)

# Step 2: Define the Quaternion RoPE rotation operator
THETA = np.pi / 16  # A fixed parameter for the encoding
# A fixed, normalized rotation axis u
U_VEC = np.array([1.0, -2.0, 3.0])
U_VEC = U_VEC / np.linalg.norm(U_VEC)
U1, U2, U3 = U_VEC

def get_R(m):
    """Generates the rotation quaternion R_m for position m."""
    angle = m * THETA
    w = np.cos(angle)
    s = np.sin(angle)
    return Quaternion(w, s * U1, s * U2, s * U3)

def analyze_statements():
    """Analyzes each statement and determines if it is true."""
    correct_statements = []
    
    # Define random quaternions for testing
    q = Quaternion(0.5, 1.2, -0.3, 0.8)
    k = Quaternion(-0.1, 0.4, 1.5, -2.1)
    
    # A) The inner product ⟨R_m(q), R_n(k)⟩ depends only on |m-n|
    m, n = 5, 8
    # The property holds if ⟨R_m(q), R_n(k)⟩ = ⟨q, R_{n-m}(k)⟩
    lhs = inner_product(get_R(m) * q, get_R(n) * k)
    rhs = inner_product(q, get_R(n-m) * k)
    if np.isclose(lhs, rhs):
        correct_statements.append('A')

    # B) R_m(q)R_n(k) = R_{m+p}(q)R_{n+p}(k) for any p
    m, n, p = 2, 3, 4
    lhs_B = (get_R(m) * q) * (get_R(n) * k)
    rhs_B = (get_R(m+p) * q) * (get_R(n+p) * k)
    if not lhs_B.is_close(rhs_B):
        # Statement is false, do nothing
        pass

    # C) The rotation axis (u_1,u_2,u_3) must be fixed for all positions
    # This is true by definition for RoPE-like relative embeddings.
    # The composition property R_m R_n = R_{m+n} depends on it.
    correct_statements.append('C')
        
    # D) Quaternion RoPE can encode 3D relative positions with a single rotation
    # The formulation uses a scalar position 'm', so it encodes 1D position.
    # Statement is false, do nothing
    pass

    # E) The magnitude |R_m(v)| equals |v| for all m
    m, v = 7, q
    lhs_E = (get_R(m) * v).norm
    rhs_E = v.norm
    if np.isclose(lhs_E, rhs_E):
        correct_statements.append('E')

    # F) R_m(αv) = αR_m(v) for scalar α
    m, v, alpha = 3, q, 2.5
    lhs_F = get_R(m) * (alpha * v)
    rhs_F = alpha * (get_R(m) * v)
    if lhs_F.is_close(rhs_F):
        correct_statements.append('F')

    # G) The quaternion rotation preserves orthogonality between vectors
    # Test with two orthogonal vectors q_o, k_o
    q_o = Quaternion(1, 1, 0, 0)
    k_o = Quaternion(0, 0, 1, -1) # <q_o, k_o> = 0
    m = 5
    # Check if <R_m(q_o), R_m(k_o)> is still 0
    inner_prod_rotated = inner_product(get_R(m) * q_o, get_R(m) * k_o)
    if np.isclose(inner_prod_rotated, 0):
        correct_statements.append('G')

    # H) The composition R_m ∘ R_n equals R_{m+n}
    m, n, v = 3, 4, q
    lhs_H = get_R(m) * (get_R(n) * v) # (R_m o R_n)(v)
    rhs_H = get_R(m+n) * v
    if lhs_H.is_close(rhs_H):
        correct_statements.append('H')
        
    # J) For any quaternion vector v, (R_m ∘ R_n)(v) - R_n(R_m(v)) is always purely imaginary
    # Since R_m and R_n commute, the result is 0, which has a real part of 0.
    m, n, v = 2, 5, q
    term1 = (get_R(m) * get_R(n)) * v
    term2 = (get_R(n) * get_R(m)) * v
    diff = term1 - term2
    if np.isclose(diff.w, 0):
        correct_statements.append('J')

    # K) The double quaternion rotation R_m(R_n(v)) preserves the real component of v
    m, n, v = 3, 4, q
    v_rotated = get_R(m+n) * v
    if not np.isclose(v_rotated.w, v.w):
        # Statement is false, do nothing
        pass
        
    # L) The trace of the 4×4 matrix representation of R_m equals 4cos(mθ)
    m = 5
    R = get_R(m)
    a, b, c, d = R.coeffs
    # Matrix for left multiplication by R
    L_R = np.array([
        [a, -b, -c, -d],
        [b,  a, -d,  c],
        [c,  d,  a, -b],
        [d, -c,  b,  a]
    ])
    trace = np.trace(L_R)
    expected_trace = 4 * np.cos(m * THETA)
    if np.isclose(trace, expected_trace):
        correct_statements.append('L')

    # M) The commutator [R_m, R_n] = 0 if and only if the rotation axis is parallel to (1,1,1)
    # The commutator is always 0 for any axis, so the "if and only if" is false.
    pass

    # N) For any unit quaternion v, lim_{m→∞} |R_m(v)| converges to 1/√2
    # From (E), |R_m(v)| = |v|. For a unit quaternion, this is 1. The limit is 1.
    pass

    # Final Step: Sort and print the results
    correct_statements.sort()
    final_answer = "".join(correct_statements)
    print("Analysis complete.")
    print(f"The true statements are: {', '.join(correct_statements)}")
    print("Final sorted list of correct statement letters:")
    print(final_answer)

if __name__ == '__main__':
    analyze_statements()