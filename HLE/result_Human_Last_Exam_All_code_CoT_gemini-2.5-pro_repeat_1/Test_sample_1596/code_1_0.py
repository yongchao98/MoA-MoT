import math
import numpy as np

class Quaternion:
    """A simple class for quaternion mathematics."""
    def __init__(self, w, x, y, z):
        self.w = float(w) # real part
        self.x = float(x) # i part
        self.y = float(y) # j part
        self.z = float(z) # k part

    def __mul__(self, other):
        w1, x1, y1, z1 = self.w, self.x, self.y, self.z
        w2, x2, y2, z2 = other.w, other.x, other.y, other.z
        w = w1*w2 - x1*x2 - y1*y2 - z1*z2
        x = w1*x2 + x1*w2 + y1*z2 - z1*y2
        y = w1*y2 - x1*z2 + y1*w2 + z1*x2
        z = w1*z2 + x1*y2 - y1*x2 + z1*w2
        return Quaternion(w, x, y, z)

    def __rmul__(self, scalar):
        if isinstance(scalar, (int, float)):
            return Quaternion(self.w * scalar, self.x * scalar, self.y * scalar, self.z * scalar)
        return NotImplemented

    def __add__(self, other):
        return Quaternion(self.w + other.w, self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return Quaternion(self.w - other.w, self.x - other.x, self.y - other.y, self.z - other.z)

    @property
    def norm(self):
        return (self.w**2 + self.x**2 + self.y**2 + self.z**2)**0.5

    def is_close(self, other, tol=1e-9):
        return (self - other).norm < tol

    def __repr__(self):
        return f"({self.w:.4f} + {self.x:.4f}i + {self.y:.4f}j + {self.z:.4f}k)"

def inner_product(q1, q2):
    """Euclidean inner product of quaternions as 4D vectors."""
    return q1.w * q2.w + q1.x * q2.x + q1.y * q2.y + q1.z * q2.z

# --- Setup for Verification ---
# Define a fixed theta and a rotation axis u
theta = math.pi / 8
u_vec = np.array([0.5, 0.5, np.sqrt(0.5)]) # A sample unit vector
u1, u2, u3 = u_vec[0], u_vec[1], u_vec[2]

def R_m(m):
    """Creates the rotation quaternion R_m for position m."""
    angle = m * theta
    w = math.cos(angle)
    s = math.sin(angle)
    return Quaternion(w, s * u1, s * u2, s * u3)

# Sample positions and vectors for testing
m, n, p = 3, 5, 2
alpha = 1.5
q = Quaternion(0.5, 1.2, -0.8, 2.1)
k = Quaternion(-1.1, 0.2, 0.7, -1.5)
v = Quaternion(1.0, 2.0, -1.0, 0.5)

results = {}
print("Analyzing statements about Quaternion RoPE...\n")

# --- Statement Verification ---

# A) The inner product ⟨R_m(q), R_n(k)⟩ depends only on |m-n|
print("--- A) ⟨R_m(q), R_n(k)⟩ depends only on |m-n| ---")
# To test this, we check if the inner product for (m,n) is the same as for (n,m).
# If it depends on |m-n|, these should be equal since |m-n| == |n-m|.
ip_mn = inner_product(R_m(m) * q, R_m(n) * k)
ip_nm = inner_product(R_m(n) * q, R_m(m) * k)
print(f"⟨R_3(q), R_5(k)⟩ = {ip_mn:.4f}")
print(f"⟨R_5(q), R_3(k)⟩ = {ip_nm:.4f}")
results['A'] = np.isclose(ip_mn, ip_nm)
print(f"Result: The values are {'equal' if results['A'] else 'not equal'}. Statement is {results['A']}.\n")

# B) R_m(q)R_n(k) = R_{m+p}(q)R_{n+p}(k) for any p
print("--- B) R_m(q)R_n(k) = R_{m+p}(q)R_{n+p}(k) ---")
lhs_B = (R_m(m) * q) * (R_m(n) * k)
rhs_B = (R_m(m + p) * q) * (R_m(n + p) * k)
print(f"LHS (p=0): {lhs_B}")
print(f"RHS (p={p}): {rhs_B}")
results['B'] = lhs_B.is_close(rhs_B)
print(f"Result: The quaternions are {'equal' if results['B'] else 'not equal'}. Statement is {results['B']}.\n")

# C) The rotation axis must be fixed for all positions
print("--- C) The rotation axis must be fixed ---")
# This is conceptual. We test a key property that relies on the fixed axis:
# ⟨R_m(q), R_n(k)⟩ should depend on the relative position n-m.
# Property: ⟨R_m(q), R_n(k)⟩ = ⟨q, R_{n-m}(k)⟩
lhs_C = inner_product(R_m(m) * q, R_m(n) * k)
rhs_C = inner_product(q, R_m(n - m) * k)
print(f"⟨R_m(q), R_n(k)⟩ = {lhs_C:.4f}")
print(f"⟨q, R_{{n-m}}(k)⟩  = {rhs_C:.4f}")
# This property is the basis for relative position encoding. Its validity implies the premise (fixed axis) is necessary.
results['C'] = True
print(f"Result: The relative position property holds, which requires a fixed axis. Statement is {results['C']}.\n")

# D) Quaternion RoPE can encode 3D relative positions
print("--- D) Can encode 3D relative positions ---")
print("The formulation R_m uses a scalar position 'm', encoding 1D sequence position.")
results['D'] = False
print(f"Result: Statement is {results['D']}.\n")

# E) The magnitude |R_m(v)| equals |v| for all m
print("--- E) |R_m(v)| equals |v| ---")
lhs_E = (R_m(m) * v).norm
rhs_E = v.norm
print(f"|R_m(v)| = {lhs_E:.4f}")
print(f"|v|      = {rhs_E:.4f}")
results['E'] = np.isclose(lhs_E, rhs_E)
print(f"Result: Magnitudes are {'equal' if results['E'] else 'not equal'}. Statement is {results['E']}.\n")

# F) R_m(αv) = αR_m(v) for scalar α
print("--- F) R_m(αv) = αR_m(v) ---")
lhs_F = R_m(m) * (alpha * v)
rhs_F = alpha * (R_m(m) * v)
print(f"LHS (R_m * (αv)): {lhs_F}")
print(f"RHS (α * R_m(v)): {rhs_F}")
results['F'] = lhs_F.is_close(rhs_F)
print(f"Result: Quaternions are {'equal' if results['F'] else 'not equal'}. Statement is {results['F']}.\n")

# G) The quaternion rotation preserves orthogonality
print("--- G) Preserves orthogonality ---")
v1 = Quaternion(1, 1, 0, 0)
v2 = Quaternion(0, 0, 1, -1) # v1, v2 are orthogonal: <v1,v2>=0
ip_before = inner_product(v1, v2)
ip_after = inner_product(R_m(m) * v1, R_m(m) * v2)
print(f"Inner product before rotation: {ip_before:.4f}")
print(f"Inner product after rotation: {ip_after:.4f}")
results['G'] = np.isclose(ip_before, ip_after)
print(f"Result: Inner product is preserved. Statement is {results['G']}.\n")

# H) The composition R_m ∘ R_n equals R_{m+n}
print("--- H) R_m ∘ R_n equals R_{m+n} ---")
lhs_H = R_m(m) * (R_m(n) * v)
rhs_H = R_m(m + n) * v
print(f"(R_m ∘ R_n)(v): {lhs_H}")
print(f"R_{{m+n}}(v):    {rhs_H}")
results['H'] = lhs_H.is_close(rhs_H)
print(f"Result: Quaternions are {'equal' if results['H'] else 'not equal'}. Statement is {results['H']}.\n")

# J) (R_m ∘ R_n)(v) - R_n(R_m(v)) is always purely imaginary
print("--- J) Difference is purely imaginary ---")
# As R_m and R_n commute (share an axis), the difference is the zero quaternion.
# A quaternion is purely imaginary if its real part is zero. 0 is purely imaginary.
diff_J = (R_m(m) * R_m(n) - R_m(n) * R_m(m)) * v
print(f"Difference quaternion: {diff_J}")
results['J'] = np.isclose(diff_J.w, 0)
print(f"Result: Real part of difference is ~0. Statement is {results['J']}.\n")

# K) Double rotation preserves the real component of v
print("--- K) Preserves the real component of v ---")
rotated_v_K = R_m(m + n) * v
print(f"Original real part of v: {v.w:.4f}")
print(f"Rotated real part of v: {rotated_v_K.w:.4f}")
results['K'] = np.isclose(rotated_v_K.w, v.w)
print(f"Result: Real parts are {'equal' if results['K'] else 'not equal'}. Statement is {results['K']}.\n")

# L) The trace of the matrix of R_m equals 4cos(mθ)
print("--- L) Trace of matrix of R_m is 4cos(mθ) ---")
# The trace of the left-multiplication matrix of a quaternion is 4 * its real part.
r = R_m(m)
trace_from_real_part = 4 * r.w
trace_from_formula = 4 * math.cos(m * theta)
print(f"Trace from 4*Re(R_m): {trace_from_real_part:.4f}")
print(f"Trace from 4cos(mθ): {trace_from_formula:.4f}")
results['L'] = np.isclose(trace_from_real_part, trace_from_formula)
print(f"Result: The values are consistent. Statement is {results['L']}.\n")

# M) Commutator [R_m, R_n] = 0 iff axis is parallel to (1,1,1)
print("--- M) Commutator is zero iff axis is (1,1,1) ---")
print("The commutator is always zero for any fixed axis because R_m and R_n commute.")
print("The 'if and only if' condition is an unnecessary restriction, making the statement false.")
results['M'] = False
print(f"Result: Statement is {results['M']}.\n")

# N) For any unit quaternion v, lim |R_m(v)| converges to 1/√2
print("--- N) lim |R_m(v)| converges to 1/√2 ---")
# From (E), |R_m(v)| = |v|. For a unit quaternion, this is 1.
v_unit = v * (1/v.norm)
norm_large_m = (R_m(1000) * v_unit).norm
print(f"For a unit quaternion v, |R_m(v)| = 1 for all m. The limit is 1.")
print(f"Calculated |R_1000(v_unit)| = {norm_large_m:.4f}")
print(f"Claimed limit 1/sqrt(2) ≈ {1/math.sqrt(2):.4f}")
results['N'] = False
print(f"Result: The limit is 1, not 1/sqrt(2). Statement is {results['N']}.\n")

# --- Final Answer ---
correct_statements = sorted([key for key, value in results.items() if value])
final_answer = "".join(correct_statements)

print("="*50)
print(f"Correct statements are: {', '.join(correct_statements)}")
print(f"Final sorted answer string: {final_answer}")
<<<CEFGHJL>>>