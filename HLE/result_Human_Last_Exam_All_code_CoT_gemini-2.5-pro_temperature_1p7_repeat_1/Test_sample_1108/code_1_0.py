import sympy as sp

# Step 1: Define symbolic variables and parameters
S, F = sp.symbols('S F')
a, b, c, h, m = 2, 1, 1, 1, 1

# Step 2: Define the system's rate functions
G = S * (h - m * S / F)
H = F * (a - b * F - c * S)

# Step 3: Define the non-trivial equilibrium point found manually
S0 = 1
F0 = 1

# Step 4: Compute the Jacobian matrix of the system [G, H] with respect to [S, F]
J = sp.Matrix([G, H]).jacobian([S, F])

# Step 5: Evaluate the Jacobian at the equilibrium point (S0, F0)
J_at_eq = J.subs({S: S0, F: F0})

# Step 6: Extract the coefficients a_ij from the evaluated Jacobian matrix
a11 = J_at_eq[0, 0]
a12 = J_at_eq[0, 1]
a21 = J_at_eq[1, 0]
a22 = J_at_eq[1, 1]

# Step 7: The coefficients b_ij are the values of the rate functions at equilibrium,
# which are zero by definition.
b11 = G.subs({S: S0, F: F0})
b22 = H.subs({S: S0, F: F0})

# Step 8: Print the final coefficients as requested.
print(f"The non-trivial equilibrium point is (S0, F0) = ({S0}, {F0}).")
print("The coefficients of the linearized system are:")
print(f"a_11 = {a11}")
print(f"a_12 = {a12}")
print(f"a_21 = {a21}")
print(f"a_22 = {a22}")
print(f"b_11 = {b11}")
print(f"b_22 = {b22}")
