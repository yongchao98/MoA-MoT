import sympy

# --- Step 1: Define the system symbolically ---
# Define symbolic variables for S, F, and the parameters
S, F = sympy.symbols('S F')
a, b, c, h, m = sympy.symbols('a b c h m')

# Define the system of ODEs
f_S = S * (h - m * S / F)
f_F = F * (a - b * F - c * S)

# Substitute the given parameter values
params = {a: 2, b: 1, c: 1, h: 1, m: 1}
f_S_subs = f_S.subs(params)
f_F_subs = f_F.subs(params)

# --- Step 2: Find the non-trivial equilibrium point ---
# From manual derivation, we know the equilibrium point is (1, 1)
S_e, F_e = 1, 1

# --- Step 3 & 4: Compute the Jacobian and evaluate it ---
# Create a matrix for the system functions
F_matrix = sympy.Matrix([f_S_subs, f_F_subs])

# Create a matrix for the variables
X_matrix = sympy.Matrix([S, F])

# Compute the Jacobian matrix
J = F_matrix.jacobian(X_matrix)

# Substitute the equilibrium point into the Jacobian
J_eq = J.subs({S: S_e, F: F_e})

# --- Step 5: Extract coefficients ---
# The elements of the Jacobian are the a_ij coefficients
a11 = J_eq[0, 0]
a12 = J_eq[0, 1]
a21 = J_eq[1, 0]
a22 = J_eq[1, 1]

# The linearization around an equilibrium point results in a homogeneous system,
# so the constant vector 'b' is zero.
b11 = 0
b22 = 0

# Print the final results
print(f"a_11 = {a11}")
print(f"a_12 = {a12}")
print(f"a_21 = {a21}")
print(f"a_22 = {a22}")
print(f"b_11 = {b11}")
print(f"b_22 = {b22}")

# Format the final answer string as requested
final_answer_str = f"a_11 = {a11}, a_12 = {a12}, a_21 = {a21}, a_22 = {a22}, b_11 = {b11}, b_22 = {b22}"
print(f"<<<{final_answer_str}>>>")