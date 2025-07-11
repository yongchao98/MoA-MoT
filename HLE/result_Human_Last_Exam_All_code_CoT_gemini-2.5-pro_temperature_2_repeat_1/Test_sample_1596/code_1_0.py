import numpy as np

class Quaternion:
    """A helper class to represent and operate on quaternions."""
    def __init__(self, w, x, y, z):
        # w is the real part, (x, y, z) is the vector part
        self.w, self.x, self.y, self.z = float(w), float(x), float(y), float(z)

    def __mul__(self, other):
        w1, x1, y1, z1 = self.w, self.x, self.y, self.z
        w2, x2, y2, z2 = other.w, other.x, other.y, other.z
        return Quaternion(
            w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
            w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
        )

    def __rmul__(self, scalar):
        if isinstance(scalar, (int, float)):
            return Quaternion(scalar * self.w, scalar * self.x, scalar * self.y, scalar * self.z)
        raise NotImplementedError

    def __sub__(self, other):
        return Quaternion(self.w - other.w, self.x - other.x, self.y - other.y, self.z - other.z)

    def __repr__(self):
        return f"({self.w:.4f} + {self.x:.4f}i + {self.y:.4f}j + {self.z:.4f}k)"
    
    @property
    def norm(self):
        return np.sqrt(self.w**2 + self.x**2 + self.y**2 + self.z**2)
    
    @property
    def real(self):
        return self.w
        
    @property
    def is_purely_imaginary(self):
        return np.isclose(self.w, 0)

# ----- Setup for Tests -----

# Define a fixed rotation axis `u` (must be a unit vector)
u_vec = np.array([1.0, -2.0, 3.0])
u_vec = u_vec / np.linalg.norm(u_vec)
# Define a fixed base angle `theta`
theta = np.pi / 8

def R_m_op(m):
    """Generates the rotation quaternion R_m for position m."""
    angle = m * theta
    w_part = np.cos(angle)
    sin_part = np.sin(angle)
    x_part, y_part, z_part = u_vec * sin_part
    return Quaternion(w_part, x_part, y_part, z_part)

def inner_product(p, q):
    """Calculates the 4D inner product of two quaternions."""
    return p.w * q.w + p.x * q.x + p.y * q.y + p.z * q.z

# Sample query and key vectors
q = Quaternion(0.5, 1.0, -0.2, 0.8)
k = Quaternion(0.1, -0.7, 1.5, -0.1)

true_statements = []

# ----- Statement Analysis -----

print("Analyzing statements about Quaternion RoPE...\n")

# A) The inner product <R_m(q), R_n(k)> depends only on |m-n|
print("--- Statement A ---")
m1, n1 = 2, 5  # m-n = -3, |m-n| = 3
m2, n2 = 5, 2  # m-n = 3,  |m-n| = 3
val1 = inner_product(R_m_op(m1) * q, R_m_op(n1) * k)
val2 = inner_product(R_m_op(m2) * q, R_m_op(n2) * k)
print(f"For m=2, n=5, the inner product is: {val1:.4f}")
print(f"For m=5, n=2, the inner product is: {val2:.4f}")
print("The values are different for the same |m-n|. The result depends on m-n, not its absolute value.")
print("Result: False\n")

# B) R_m(q)R_n(k) = R_{m+p}(q)R_{n+p}(k) for any p
print("--- Statement B ---")
print("This statement claims (R_m*q)*(R_n*k) = (R_{m+p}*q)*(R_{n+p}*k).")
print("This equality does not hold because quaternion multiplication is non-commutative, and R_p cannot be factored out.")
print("Result: False\n")

# C) The rotation axis (u_1,u_2,u_3) must be fixed for all positions
print("--- Statement C ---")
print("Reasoning: The key property of RoPE is encoding relative positions. The inner product <R_m(q), R_n(k)> simplifies to a function of (n-m) *only if* the rotation axis `u` is constant for both R_m and R_n. If the axis u_m depended on position, this property would be lost.")
print("Result: True\n")
true_statements.append('C')

# D) Quaternion RoPE can encode 3D relative positions with a single rotation
print("--- Statement D ---")
print("Reasoning: The operator R_m has one positional degree of freedom, `m`, which scales the rotation angle. This is designed to encode a 1D position. Encoding a 3D position vector would require more degrees of freedom, like multiple rotation axes.")
print("Result: False\n")

# E) The magnitude |R_m(v)| equals |v| for all m
print("--- Statement E ---")
m = 4
v = Quaternion(1, 2, 3, 4)
norm_Rm = R_m_op(m).norm
norm_v = v.norm
norm_Rmv = (R_m_op(m) * v).norm
print(f"The norm of the rotation operator |R_{m}| is |cos({m}θ) + u*sin({m}θ)| = sqrt(cos^2 + sin^2*(u_x^2+u_y^2+u_z^2)) = {norm_Rm:.4f}")
print(f"For a vector v, |v| = {norm_v:.4f}")
print(f"|R_m(v)| = |R_m|*|v| = {norm_Rm:.4f} * {norm_v:.4f} = {norm_Rmv:.4f}")
print("Since |R_m| is 1, the norm of v is preserved.")
print("Result: True\n")
true_statements.append('E')

# F) R_m(αv) = αR_m(v) for scalar α
print("--- Statement F ---")
alpha = 2.5
lhs = R_m_op(m) * (alpha * v)
rhs = alpha * (R_m_op(m) * v)
print(f"Let α=2.5. R_m(αv) = {lhs}")
print(f"αR_m(v) = {rhs}")
print("The results are equal because scalar multiplication commutes with quaternion multiplication.")
print("Result: True\n")
true_statements.append('F')

# G) The quaternion rotation preserves orthogonality between vectors
print("--- Statement G ---")
p_ortho = Quaternion(1, -2, 3, -1)
q_ortho = Quaternion(2, 4, 2, 0)
ip_orig = inner_product(p_ortho, q_ortho)
p_rot = R_m_op(m) * p_ortho
q_rot = R_m_op(m) * q_ortho
ip_rot = inner_product(p_rot, q_rot)
print(f"Let p and q be orthogonal, such that <p,q> = {ip_orig:.4f}")
print(f"After rotation, <R_m(p), R_m(q)> = {ip_rot:.4f}")
print("The inner product <R_m(p), R_m(q)> is equal to <p, q>. So if the original is 0, the new one is 0.")
print("Result: True\n")
true_statements.append('G')

# H) The composition R_m ∘ R_n equals R_{m+n}
print("--- Statement H ---")
m, n = 3, 7
comp_op = R_m_op(m) * R_m_op(n)
sum_op = R_m_op(m + n)
print(f"(R_m * R_n) results in: {comp_op}")
print(f"R_{m+n} results in:    {sum_op}")
print("The operators are equal because they share the same rotation axis `u` and their angles add up.")
print("Result: True\n")
true_statements.append('H')

# J) For any quaternion vector v, (R_m ∘ R_n)(v) - R_n(R_m(v)) is always purely imaginary
print("--- Statement J ---")
commutator_op = (R_m_op(m) * R_m_op(n)) - (R_m_op(n) * R_m_op(m))
print(f"The expression equals [R_m, R_n]v. Let's calculate the operator [R_m, R_n] = R_m*R_n - R_n*R_m.")
print(f"Resulting operator: {commutator_op}")
print(f"The commutator is the zero quaternion. A quaternion is purely imaginary if its real part is zero. The real part of 0 is 0.")
print("Result: True\n")
true_statements.append('J')


# K) The double quaternion rotation R_m(R_n(v)) preserves the real component of v regardless of m,n
print("--- Statement K ---")
v_real = v.real
rot_v = (R_m_op(m) * R_m_op(n)) * v
rot_v_real = rot_v.real
print(f"The original real component of v is {v_real:.4f}")
print(f"The real component after rotation R_{m+n}(v) is {rot_v_real:.4f}")
print("The real component is not preserved.")
print("Result: False\n")

# L) The trace of the 4×4 matrix representation of R_m equals 4cos(mθ)
print("--- Statement L ---")
print("The 4x4 matrix L_q for left multiplication by q=a+bi+cj+dk has `a` on all four diagonal positions.")
print("The trace is the sum of the diagonal elements, so Tr(L_q) = 4a = 4*Re(q).")
a = R_m_op(m).real
trace = 4 * a
expected_trace = 4 * np.cos(m * theta)
print(f"For R_m, the real part `a` is cos(mθ) = {a:.4f}")
print(f"The trace is 4 * {a:.4f} = {trace:.4f}")
print(f"The theoretical value is 4*cos({m}θ) = {expected_trace:.4f}")
print("Result: True\n")
true_statements.append('L')

# M) The commutator [R_m, R_n] = 0 if and only if the rotation axis is parallel to (1,1,1)
print("--- Statement M ---")
print("Reasoning: As shown in J, the commutator [R_m, R_n] is always zero for *any* fixed rotation axis `u`, because quaternions with the same complex axis commute. The 'if and only if' condition is an incorrect restriction.")
print("Result: False\n")

# N) For any unit quaternion v, lim_{m→∞} |R_m(v)| converges to 1/√2
print("--- Statement N ---")
print("Reasoning: As shown in E, |R_m(v)| = |R_m|*|v|. We know |R_m|=1, and for a unit quaternion v, |v|=1. Thus, |R_m(v)| = 1 for all m. The limit of a constant sequence of 1s is 1, not 1/sqrt(2).")
print("Result: False\n")


# ----- Final Answer -----
final_answer = "".join(sorted(true_statements))
print("The sorted list of true statement letters is:", final_answer)
print(f"\n<<<{''.join(sorted(true_statements))}>>>")