import sympy

# Define symbols for the variables and parameters
S, F, a, b, c, h, m = sympy.symbols('S F a b c h m')

# Define the system of equations
dS_dt = S * (h - m * S / F)
dF_dt = F * (a - b*F - c*S)

# Substitute the given parameter values
params = {a: 2, b: 1, c: 1, h: 1, m: 1}
dS_dt_sub = dS_dt.subs(params)
dF_dt_sub = dF_dt.subs(params)

# Solve for the non-trivial equilibrium point (S > 0, F > 0)
# From dS/dt = 0, we have h - m*S/F = 0 => S = F * h / m
# Substitute into dF/dt = 0: a - b*F - c*(F*h/m) = 0
# Solve for F: a = F * (b + c*h/m) => F = a / (b + c*h/m)
a_val, b_val, c_val, h_val, m_val = params.values()
F_e = a_val / (b_val + c_val * h_val / m_val)
S_e = F_e * h_val / m_val

print(f"The non-trivial equilibrium solution is (S_e, F_e) = ({S_e}, {F_e})\n")

# Compute the Jacobian matrix
J = sympy.Matrix([dS_dt_sub, dF_dt_sub]).jacobian([S, F])
# print("Jacobian Matrix J(S, F):")
# sympy.pprint(J)

# Evaluate the Jacobian at the equilibrium point
J_eq = J.subs({S: S_e, F: F_e})
# print(f"\nJacobian at equilibrium point ({S_e}, {F_e}):")
# sympy.pprint(J_eq)

a11 = J_eq[0, 0]
a12 = J_eq[0, 1]
a21 = J_eq[1, 0]
a22 = J_eq[1, 1]

# The constant terms b11 and b22 are zero for linearization at an equilibrium point.
b11 = 0
b22 = 0

# Print the final coefficients
print("The coefficients of the linearized system are:")
print(f"a_11 = {a11}")
print(f"a_12 = {a12}")
print(f"a_21 = {a21}")
print(f"a_22 = {a22}")
print(f"b_11 = {b11}")
print(f"b_22 = {b22}")
