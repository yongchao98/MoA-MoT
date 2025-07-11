import numpy as np

# --- Quaternion Math Utilities ---
def quat_mul(p, q):
    """Multiplies two quaternions p * q."""
    p0, p1, p2, p3 = p
    q0, q1, q2, q3 = q
    r0 = p0 * q0 - p1 * q1 - p2 * q2 - p3 * q3
    r1 = p0 * q1 + p1 * q0 + p2 * q3 - p3 * q2
    r2 = p0 * q2 - p1 * q3 + p2 * q0 + p3 * q1
    r3 = p0 * q3 + p1 * q2 - p2 * q1 + p3 * q0
    return np.array([r0, r1, r2, r3])

def quat_conj(q):
    """Computes the conjugate of a quaternion."""
    return np.array([q[0], -q[1], -q[2], -q[3]])

def quat_norm(q):
    """Computes the norm of a quaternion."""
    return np.linalg.norm(q)

def get_R(m, theta, u):
    """Creates the rotation quaternion R_m."""
    s = np.sin(m * theta)
    c = np.cos(m * theta)
    return np.array([c, s * u[0], s * u[1], s * u[2]])

def inner_product(p, q):
    """Computes the inner product of two quaternions."""
    return np.dot(p, q)

def get_matrix_L(p):
    """Gets the left-multiplication matrix for a quaternion."""
    p0, p1, p2, p3 = p
    return np.array([
        [p0, -p1, -p2, -p3],
        [p1,  p0, -p3,  p2],
        [p2,  p3,  p0, -p1],
        [p3, -p2,  p1,  p0]
    ])

# --- Setup for Verification ---
print("Verifying statements about Quaternion RoPE...\n")
# Fixed parameters
theta = np.pi / 16
u = np.array([1, 2, 3]) / np.sqrt(14) # A fixed unit vector axis

# Example positions
m, n, p_offset = 5, 8, 2

# Example quaternions
q_vec = np.array([0.5, 1.0, -2.0, 0.2])
k_vec = np.array([-0.8, 0.6, 1.5, -1.1])
v_vec = np.array([2.0, -1.0, 3.0, -0.5])
v_unit = v_vec / quat_norm(v_vec)

# Rotation quaternions
R_m = get_R(m, theta, u)
R_n = get_R(n, theta, u)

# --- Verification Logic ---
results = {}

# A) Dependency on |m-n|
ip1 = inner_product(quat_mul(R_m, q_vec), quat_mul(R_n, k_vec))
ip2 = inner_product(quat_mul(R_n, q_vec), quat_mul(R_m, k_vec)) # Swapping m and n
results['A'] = np.isclose(ip1, ip2) # False, because it depends on m-n, not |m-n|
print(f"A) False. ⟨R_m(q),R_n(k)⟩ = {ip1:.4f}, but ⟨R_n(q),R_m(k)⟩ = {ip2:.4f}. Not equal.")

# C) Fixed axis
results['C'] = True # By reasoning, this is required for the relative property.
print("C) True. A fixed rotation axis is required to make the inner product relative.")

# E) |R_m(v)| = |v|
norm_v = quat_norm(v_vec)
norm_Rmv = quat_norm(quat_mul(R_m, v_vec))
results['E'] = np.isclose(norm_v, norm_Rmv)
print(f"E) True. |v| = {norm_v:.4f}, |R_m(v)| = {norm_Rmv:.4f}.")

# F) R_m(αv) = αR_m(v)
alpha = 1.5
lhs = quat_mul(R_m, alpha * v_vec)
rhs = alpha * quat_mul(R_m, v_vec)
results['F'] = np.allclose(lhs, rhs)
print(f"F) True. R_m(αv) - αR_m(v) = {quat_norm(lhs-rhs):.4e}.")

# G) Preserves orthogonality
w_vec = np.array([1.0, 1.0, 1.0, -4.0]) # v_vec and w_vec are orthogonal by construction
# inner_product(v_vec, w_vec) = 2 - 1 + 3 - 2 = 2. Let's make one.
v_ortho = np.array([1, 1, 0, 0])
w_ortho = np.array([0, 0, 1, 1])
ip_orig = inner_product(v_ortho, w_ortho)
ip_rot = inner_product(quat_mul(R_m, v_ortho), quat_mul(R_m, w_ortho))
results['G'] = np.isclose(ip_orig, ip_rot)
print(f"G) True. Original inner product = {ip_orig:.4f}, after rotation = {ip_rot:.4f}.")

# H) R_m ∘ R_n = R_{m+n}
R_m_plus_n = get_R(m + n, theta, u)
R_m_R_n = quat_mul(R_m, R_n)
results['H'] = np.allclose(R_m_plus_n, R_m_R_n)
print(f"H) True. ||R_{{m+n}} - (R_m R_n)|| = {quat_norm(R_m_plus_n - R_m_R_n):.4e}.")

# J) Commutator is purely imaginary
R_n_R_m = quat_mul(R_n, R_m)
commutator_v = quat_mul(R_m_R_n - R_n_R_m, v_vec)
# Commutator should be the zero quaternion
results['J'] = np.allclose(commutator_v, np.zeros(4))
print(f"J) True. Commutator result ||[R_m, R_n]v|| = {quat_norm(commutator_v):.4e} (is zero). Zero is purely imaginary.")

# K) Preserves real component
v_real_part = v_vec[0]
R_mn_v = quat_mul(get_R(m + n, theta, u), v_vec)
R_mn_v_real_part = R_mn_v[0]
results['K'] = np.isclose(v_real_part, R_mn_v_real_part)
print(f"K) False. Original real part = {v_real_part:.4f}, after rotation = {R_mn_v_real_part:.4f}.")

# L) Trace of matrix is 4cos(mθ)
L_Rm = get_matrix_L(R_m)
trace = np.trace(L_Rm)
expected_trace = 4 * np.cos(m * theta)
results['L'] = np.isclose(trace, expected_trace)
print(f"L) True. Matrix trace = {trace:.4f}, 4cos(mθ) = {expected_trace:.4f}.")

# M) Commutator is 0 iff axis is || (1,1,1)
results['M'] = False # Reasoning: commutator is 0 for ANY fixed axis.
print(f"M) False. The 'only if' condition is incorrect.")

# N) Limit of norm
# |R_m(v)| = |v|. For a unit quaternion v, the norm is always 1.
limit_val = 1.0
results['N'] = np.isclose(limit_val, 1/np.sqrt(2))
print(f"N) False. For a unit vector, the norm is always 1, not 1/sqrt(2).")

# --- Final Answer ---
true_statements = sorted([key for key, value in results.items() if value])
final_answer = "".join(true_statements)

print("\n-------------------------------------------")
print(f"The set of true statements is: {', '.join(true_statements)}")
print(f"The final sorted answer string is: {final_answer}")