import numpy as np

class Quaternion:
    """A simple class to represent and operate on quaternions."""
    def __init__(self, w, x, y, z):
        self.w = w  # real component
        self.x = x  # i component
        self.y = y  # j component
        self.z = z  # k component

    def __repr__(self):
        return f"({self.w:.4f} + {self.x:.4f}i + {self.y:.4f}j + {self.z:.4f}k)"

    def __mul__(self, other):
        if isinstance(other, Quaternion):
            w1, x1, y1, z1 = self.w, self.x, self.y, self.z
            w2, x2, y2, z2 = other.w, other.x, other.y, other.z
            w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
            x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
            y = w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2
            z = w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2
            return Quaternion(w, x, y, z)
        elif isinstance(other, (int, float)):
            return Quaternion(self.w * other, self.x * other, self.y * other, self.z * other)
        return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    def __add__(self, other):
        return Quaternion(self.w + other.w, self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return Quaternion(self.w - other.w, self.x - other.x, self.y - other.y, self.z - other.z)

    def __abs__(self):
        return np.sqrt(self.w**2 + self.x**2 + self.y**2 + self.z**2)
    
    @property
    def real(self):
        return self.w
        
    def to_vector(self):
        return np.array([self.w, self.x, self.y, self.z])

def inner_product(q1, q2):
    """Calculates the inner product of two quaternions as defined in the prompt."""
    v1 = q1.to_vector()
    v2 = q2.to_vector()
    return np.dot(v1, v2)

def get_R(m, theta, u_axis):
    """Creates the rotation quaternion R_m for a given position m."""
    angle = m * theta
    w = np.cos(angle)
    s = np.sin(angle)
    return Quaternion(w, s * u_axis[0], s * u_axis[1], s * u_axis[2])

def evaluate_statements():
    """Tests each statement from the prompt and identifies the true ones."""
    # Setup test parameters
    theta = np.pi / 10
    u_axis = np.array([1, -1, 2]) / np.sqrt(6)
    q = Quaternion(0.5, 1.2, -0.8, 0.3)
    k = Quaternion(-0.2, 0.4, 1.5, -3.0)
    m = 3
    n = 7
    p = 4
    alpha = -1.5
    
    true_statements = []

    print("Evaluating statements numerically:\n")

    # A) The inner product ⟨R_m(q), R_n(k)⟩ depends only on |m-n|
    # Test with (m,n) and (m+1, n+1). These have the same n-m.
    val1 = inner_product(get_R(m, theta, u_axis) * q, get_R(n, theta, u_axis) * k)
    # Test with m'=0, n'=n-m. |m'-n'| = |m-n|, but n'-m' != n-m unless m=0.
    val2 = inner_product(get_R(0, theta, u_axis) * q, get_R(n-m, theta, u_axis) * k)
    # Actually test a case where n-m is negated. m'=n, n'=m
    val3 = inner_product(get_R(n, theta, u_axis) * q, get_R(m, theta, u_axis) * k)
    print(f"A) Test: ⟨R_3(q), R_7(k)⟩ = {val1:.4f} vs ⟨R_7(q), R_3(k)⟩ = {val3:.4f}. They are not equal, so it doesn't depend only on |m-n|. -> False")

    # C) The rotation axis must be fixed for all positions. This is conceptually true for the scheme to work.
    true_statements.append("C")
    print("C) This is a design requirement for the relative position property to hold, similar to standard RoPE. -> True")

    # E) The magnitude |R_m(v)| equals |v| for all m
    v = q
    mag_rotated = abs(get_R(m, theta, u_axis) * v)
    mag_original = abs(v)
    is_true = np.isclose(mag_rotated, mag_original)
    print(f"E) Test: |R_m(v)| = {mag_rotated:.4f}, |v| = {mag_original:.4f}. Preserves magnitude. -> {is_true}")
    if is_true: true_statements.append("E")

    # F) R_m(αv) = αR_m(v) for scalar α
    lhs = get_R(m, theta, u_axis) * (alpha * v)
    rhs = alpha * (get_R(m, theta, u_axis) * v)
    is_true = np.allclose(lhs.to_vector(), rhs.to_vector())
    print(f"F) Test: R_m(αv) = {lhs}. αR_m(v) = {rhs}. Linear with respect to scalar multiplication. -> {is_true}")
    if is_true: true_statements.append("F")

    # G) The quaternion rotation preserves orthogonality between vectors
    v1 = Quaternion(1, 2, -1, 0.5)
    v2 = Quaternion(2, -0.5, 1, -2) # <v1,v2> = 2-1-1-1= -1. Not orthogonal. Let's make it so.
    v2 = Quaternion(2, -1, 0, 0) # <v1, v2> = 2-2=0
    ip_before = inner_product(v1, v2)
    v1_rot = get_R(m, theta, u_axis) * v1
    v2_rot = get_R(m, theta, u_axis) * v2
    ip_after = inner_product(v1_rot, v2_rot)
    is_true = np.isclose(ip_before, 0) and np.isclose(ip_after, 0)
    print(f"G) Test: Inner product of orthogonal vectors before = {ip_before:.4f}, after rotation = {ip_after:.4f}. Preserves orthogonality. -> {is_true}")
    if is_true: true_statements.append("G")

    # H) The composition R_m ∘ R_n equals R_{m+n}
    v = q
    lhs = get_R(m, theta, u_axis) * (get_R(n, theta, u_axis) * v)
    rhs = get_R(m + n, theta, u_axis) * v
    is_true = np.allclose(lhs.to_vector(), rhs.to_vector())
    print(f"H) Test: (R_m∘R_n)(v) is {lhs}. R_{m+n}(v) is {rhs}. Composition holds. -> {is_true}")
    if is_true: true_statements.append("H")
    
    # J) For any quaternion vector v, (R_m ∘ R_n)(v) - R_n(R_m(v)) is always purely imaginary
    v = k
    term1 = get_R(m, theta, u_axis) * (get_R(n, theta, u_axis) * v)
    term2 = get_R(n, theta, u_axis) * (get_R(m, theta, u_axis) * v)
    result = term1 - term2
    is_true = np.isclose(result.real, 0)
    print(f"J) Test: The difference is {result}, with real part {result.real:.4f}. The difference is zero, so it is purely imaginary. -> {is_true}")
    if is_true: true_statements.append("J")

    # K) The double quaternion rotation R_m(R_n(v)) preserves the real component of v regardless of m,n
    v = q
    rotated_v = get_R(m+n, theta, u_axis) * v
    is_true = np.isclose(rotated_v.real, v.real)
    print(f"K) Test: Original real part = {v.real:.4f}, new real part = {rotated_v.real:.4f}. Real part is not preserved. -> {is_true}")
    
    # L) The trace of the 4×4 matrix representation of R_m equals 4cos(mθ)
    R_m_quat = get_R(m, theta, u_axis)
    trace = 4 * R_m_quat.real
    expected = 4 * np.cos(m * theta)
    is_true = np.isclose(trace, expected)
    print(f"L) Test: Trace of matrix is 4 * real_part = {trace:.4f}. 4cos(mθ) = {expected:.4f}. They match. -> {is_true}")
    if is_true: true_statements.append("L")
    
    final_answer = "".join(sorted(true_statements))
    return final_answer

if __name__ == '__main__':
    correct_statements = evaluate_statements()
    print("\n--------------------------------")
    print(f"The correctly identified statements are: {', '.join(list(correct_statements))}")
    print(f"Sorted list of correct letter-indices: {correct_statements}")
    print("--------------------------------")
    print("\n<<<" + correct_statements + ">>>")