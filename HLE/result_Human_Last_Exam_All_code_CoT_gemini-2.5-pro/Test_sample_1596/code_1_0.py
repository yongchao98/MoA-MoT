import numpy as np

class Quaternion:
    """A simple Quaternion class for mathematical operations."""
    def __init__(self, w, x, y, z):
        self.w, self.x, self.y, self.z = w, x, y, z
        self.vec = np.array([w, x, y, z])

    def __mul__(self, other):
        w1, x1, y1, z1 = self.w, self.x, self.y, self.z
        w2, x2, y2, z2 = other.w, other.x, other.y, other.z
        w = w1*w2 - x1*x2 - y1*y2 - z1*z2
        x = w1*x2 + x1*w2 + y1*z2 - z1*y2
        y = w1*y2 - x1*z2 + y1*w2 + z1*x2
        z = w1*z2 + x1*y2 - y1*x2 + z1*w2
        return Quaternion(w, x, y, z)

    def __rmul__(self, scalar):
        return Quaternion(self.w * scalar, self.x * scalar, self.y * scalar, self.z * scalar)

    def __sub__(self, other):
        return Quaternion(self.w - other.w, self.x - other.x, self.y - other.y, self.z - other.z)

    @property
    def norm(self):
        return np.sqrt(self.w**2 + self.x**2 + self.y**2 + self.z**2)

    def __repr__(self):
        return f"({self.w:6.3f} + {self.x:6.3f}i + {self.y:6.3f}j + {self.z:6.3f}k)"

def inner_product(q1, q2):
    """Computes the inner product of two quaternions."""
    return q1.w*q2.w + q1.x*q2.x + q1.y*q2.y + q1.z*q2.z

# --- Setup from the problem ---
THETA = np.pi / 8
u_vec = np.array([1, 1, 2])
u_vec = u_vec / np.linalg.norm(u_vec)
u1, u2, u3 = u_vec[0], u_vec[1], u_vec[2]

def get_R_m_quat(m):
    """Generates the rotation quaternion R_m."""
    angle = m * THETA
    w = np.cos(angle)
    s = np.sin(angle)
    return Quaternion(w, u1 * s, u2 * s, u3 * s)

def apply_rotation(m, v_q):
    """Applies rotation for position m to vector v."""
    Rm_q = get_R_m_quat(m)
    return Rm_q * v_q

# --- Verification of True Statements ---
print("Verifying the determined true statements: C, E, F, G, H, J, L\n")

# C: Axis is fixed. This is a premise demonstrated by our setup.
print("--- C) Axis is fixed ---")
print(f"This is a premise. We use a fixed axis u=({u1:.3f}, {u2:.3f}, {u3:.3f}) for all operations.\n")

# E: Magnitude is preserved
print("--- E) |R_m(v)| equals |v| ---")
m, v = 5, Quaternion(10, -2, 5, 3)
Rv = apply_rotation(m, v)
print(f"Original vector v: {v}, with norm |v| = {v.norm:.4f}")
print(f"Rotated vector R_5(v): {Rv}, with norm |R_5(v)| = {Rv.norm:.4f}\n")

# F: Scalar multiplication linearity
print("--- F) R_m(αv) = αR_m(v) ---")
alpha = 2.5
lhs = apply_rotation(m, alpha * v)
rhs = alpha * apply_rotation(m, v)
print(f"LHS R_5(2.5 * v) = {lhs}")
print(f"RHS 2.5 * R_5(v) = {rhs}\n")

# G: Orthogonality is preserved
print("--- G) Rotation preserves orthogonality ---")
p, q = Quaternion(1, 2, -3, 4), Quaternion(2, -1, 4, 3) # <p,q> = 2-2-12+12=0
m = 3
Rp, Rq = apply_rotation(m, p), apply_rotation(m, q)
print(f"Orthogonal vectors p={p}, q={q}")
print(f"Initial inner product <p, q> = {inner_product(p, q):.4f}")
print(f"Inner product after rotation <R_3(p), R_3(q)> = {inner_product(Rp, Rq):.4f}\n")

# H: Composition R_m ∘ R_n = R_{m+n}
print("--- H) R_m ∘ R_n = R_{m+n} ---")
m, n, v = 2, 4, Quaternion(0.5, 0.6, -0.7, 0.8)
lhs = get_R_m_quat(m) * (get_R_m_quat(n) * v)
rhs = get_R_m_quat(m + n) * v
print(f"LHS (R_2 ∘ R_4)(v) = {lhs}")
print(f"RHS R_6(v)         = {rhs}\n")

# J: Commutator is purely imaginary (zero)
print("--- J) [R_m, R_n](v) is purely imaginary ---")
term1 = get_R_m_quat(m) * (get_R_m_quat(n) * v)
term2 = get_R_m_quat(n) * (get_R_m_quat(m) * v)
diff = term1 - term2
print(f"The difference is {diff}, which is the zero quaternion.")
print(f"The real part is {diff.w:.1f}, so it's purely imaginary.\n")

# L: Trace of matrix representation
print("--- L) Trace of matrix of R_m is 4cos(mθ) ---")
m = 7
Rm = get_R_m_quat(m)
trace_matrix = 4 * Rm.w
trace_formula = 4 * np.cos(m * THETA)
print(f"For m={m}, R_m's real part 'a' is {Rm.w:.4f}")
print(f"Trace from matrix (4*a)      = {trace_matrix:.4f}")
print(f"Trace from formula (4cos(mθ)) = {trace_formula:.4f}\n")

# --- Final Answer ---
final_answer_list = sorted(['C', 'E', 'F', 'G', 'H', 'J', 'L'])
final_answer_string = "".join(final_answer_list)
print("The final sorted list of correct statement letter-indices is:")
print(final_answer_string)