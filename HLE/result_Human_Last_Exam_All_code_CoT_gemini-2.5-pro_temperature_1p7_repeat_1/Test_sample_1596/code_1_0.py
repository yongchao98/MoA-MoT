import numpy as np

class Quaternion:
    """A class to represent and operate on quaternions."""
    def __init__(self, a, b, c, d):
        self.a = a  # real part
        self.b = b  # i part
        self.c = c  # j part
        self.d = d  # k part

    def __repr__(self):
        return f"({self.a:.2f} + {self.b:.2f}i + {self.c:.2f}j + {self.d:.2f}k)"

    def __add__(self, other):
        return Quaternion(self.a + other.a, self.b + other.b, self.c + other.c, self.d + other.d)

    def __sub__(self, other):
        return Quaternion(self.a - other.a, self.b - other.b, self.c - other.c, self.d - other.d)
        
    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return Quaternion(self.a * other, self.b * other, self.c * other, self.d * other)
        
        # Hamilton product for two quaternions
        a1, b1, c1, d1 = self.a, self.b, self.c, self.d
        a2, b2, c2, d2 = other.a, other.b, other.c, other.d
        
        new_a = a1*a2 - b1*b2 - c1*c2 - d1*d2
        new_b = a1*b2 + b1*a2 + c1*d2 - d1*c2
        new_c = a1*c2 - b1*d2 + c1*a2 + d1*b2
        new_d = a1*d2 + b1*c2 - c1*b2 + d1*a2
        return Quaternion(new_a, new_b, new_c, new_d)

    def norm(self):
        """Calculates the quaternion norm."""
        return np.sqrt(self.a**2 + self.b**2 + self.c**2 + self.d**2)

    def is_close(self, other, tol=1e-9):
        """Checks if two quaternions are approximately equal."""
        return (np.isclose(self.a, other.a, atol=tol) and
                np.isclose(self.b, other.b, atol=tol) and
                np.isclose(self.c, other.c, atol=tol) and
                np.isclose(self.d, other.d, atol=tol))

    def inner(self, other):
        """Calculates the 4D vector inner product."""
        return self.a * other.a + self.b * other.b + self.c * other.c + self.d * other.d

# ---- Verification Setup ----

# Fixed rotation parameter theta
THETA = np.pi / 8
# Fixed rotation axis u_hat (as a pure quaternion, e.g., axis [1,1,1]/sqrt(3))
AXIS_VEC = np.array([1, 1, 1]) / np.sqrt(3)
U_HAT = Quaternion(0, AXIS_VEC[0], AXIS_VEC[1], AXIS_VEC[2])

# Sample quaternion vectors for testing
Q_VEC = Quaternion(1, 2, 3, 4)
K_VEC = Quaternion(5, -6, 7, -8)
V_VEC = Quaternion(0.5, -0.2, 0.8, -0.1)

def get_R(m, theta, u):
    """Generates the rotation quaternion R_m."""
    angle = m * theta
    scalar_part = np.cos(angle)
    # Note: Quaternion multiplication with scalar is defined
    vector_part = u * np.sin(angle)
    # The sum of a scalar (promoted to quaternion) and a pure quaternion
    return Quaternion(scalar_part, 0, 0, 0) + vector_part

# ---- Statement Verification Functions ----

def check_A():
    """A) The inner product ⟨R_m(q), R_n(k)⟩ depends only on |m-n|"""
    # Analysis: The inner product depends on R_{n-m} = cos((n-m)θ) + ûsin((n-m)θ).
    # Since sin(-x) = -sin(x), the result depends on the sign of (n-m), not just its absolute value.
    return False

def check_B():
    """B) R_m(q)R_n(k) = R_{m+p}(q)R_{n+p}(k) for any p"""
    # This equality does not hold due to the nature of quaternion multiplication.
    # We can easily disprove with a numerical example.
    m, n, p = 1, 2, 3
    lhs = (get_R(m, THETA, U_HAT) * Q_VEC) * (get_R(n, THETA, U_HAT) * K_VEC)
    rhs = (get_R(m + p, THETA, U_HAT) * Q_VEC) * (get_R(n + p, THETA, U_HAT) * K_VEC)
    return lhs.is_close(rhs)

def check_C():
    """C) The rotation axis (u_1,u_2,u_3) must be fixed for all positions"""
    # The desirable properties like relative positioning (A, although it failed slightly)
    # and composition (H) rely on R_m and R_n sharing the same rotation axis.
    # This is a core assumption of the proposed RoPE-like model.
    return True

def check_D():
    """D) Quaternion RoPE can encode 3D relative positions with a single rotation"""
    # The rotation operator R_m depends only on one scalar parameter 'm'.
    # To encode a 3D position (dx, dy, dz), one would need three degrees of freedom.
    return False

def check_E():
    """E) The magnitude |R_m(v)| equals |v| for all m"""
    # |R_m * v| = |R_m| * |v|. We show |R_m| = 1.
    m = 5
    R_m = get_R(m, THETA, U_HAT)
    rotated_v = R_m * V_VEC
    return np.isclose(rotated_v.norm(), V_VEC.norm())

def check_F():
    """F) R_m(αv) = αR_m(v) for scalar α"""
    # This is a property of quaternion algebra: R * (αv) = α * (R * v).
    m, alpha = 3, 3.5
    R_m = get_R(m, THETA, U_HAT)
    lhs = R_m * (V_VEC * alpha)
    rhs = (R_m * V_VEC) * alpha
    return lhs.is_close(rhs)

def check_G():
    """G) The quaternion rotation preserves orthogonality between vectors"""
    # Left multiplication by a unit quaternion is an orthogonal transformation in R^4,
    # which preserves inner products, and thus orthogonality.
    p1 = Quaternion(1, 2, -3, 0.5)
    # Construct an orthogonal vector p2 where inner product is 0.
    # 1*a + 2*b - 3*c + 0.5*d = 0. Choose a=3, c=1.
    p2 = Quaternion(3, 0, 1, 0)
    R_m = get_R(2, THETA, U_HAT)
    # The inner product should be preserved after rotation.
    inner_before = p1.inner(p2)
    inner_after = (R_m * p1).inner(R_m * p2)
    return np.isclose(inner_before, inner_after)

def check_H():
    """H) The composition R_m ∘ R_n equals R_{m+n}"""
    # (R_m * R_n) * v = R_{m+n} * v. We check if R_m * R_n == R_{m+n}.
    m, n = 2, 3
    R_mn_compose = get_R(m, THETA, U_HAT) * get_R(n, THETA, U_HAT)
    R_m_plus_n = get_R(m + n, THETA, U_HAT)
    return R_mn_compose.is_close(R_m_plus_n)

def check_J():
    """J) For any quaternion vector v, (R_m ∘ R_n)(v) - R_n(R_m(v)) is always purely imaginary"""
    # From H, R_m and R_n commute (R_m*R_n = R_n*R_m), so the difference is the zero quaternion.
    # The zero quaternion (0 + 0i + 0j + 0k) has a real part of 0, so it is purely imaginary.
    return True

def check_K():
    """K) The double quaternion rotation R_m(R_n(v)) preserves the real component of v"""
    # Left multiplication by a general rotation quaternion mixes real and imaginary parts.
    # Real part of R*v is (a1*a2 - b1*b2 - c1*c2 - d1*d2), which is not generally a2.
    R_m = get_R(1, THETA, U_HAT)
    rotated_v = R_m * V_VEC
    return np.isclose(rotated_v.a, V_VEC.a) # Will be False

def check_L():
    """L) The trace of the 4×4 matrix representation of R_m equals 4cos(mθ)"""
    # The matrix for left multiplication has the real part 'a' on the diagonal four times.
    # The real part of R_m is cos(mθ). So the trace is 4*cos(mθ).
    return True

def check_M():
    """M) The commutator [R_m, R_n] = 0 if and only if the rotation axis is parallel to (1,1,1)"""
    # The commutator is always zero for any axis, because R_m and R_n rotate about the same axis.
    # The "if and only if" condition makes the statement false.
    return False

def check_N():
    """N) For any unit quaternion v, lim_{m→∞} |R_m(v)| converges to 1/√2"""
    # From E, |R_m(v)| = |v|. For a unit quaternion, |v|=1.
    # The magnitude is always 1, so the limit is 1.
    return False

def main():
    """Main function to check all statements and print the result."""
    all_statements = {
        'A': check_A, 'B': check_B, 'C': check_C, 'D': check_D,
        'E': check_E, 'F': check_F, 'G': check_G, 'H': check_H,
        'J': check_J, 'K': check_K, 'L': check_L, 'M': check_M,
        'N': check_N,
    }
    
    correct_statements = []
    for letter, check_func in sorted(all_statements.items()):
        if check_func():
            correct_statements.append(letter)
            
    result = "".join(correct_statements)
    print("Final sorted list of correct statement letter-indices:")
    print(result)
    print(f"<<<{result}>>>")

if __name__ == "__main__":
    main()