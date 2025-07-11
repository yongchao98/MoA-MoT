import sympy

# This script finds the linearization of the given predator-prey system at its non-trivial equilibrium.

# 1. Define symbolic variables and the system equations
S, F = sympy.symbols('S F')
a, b, c, h, m = sympy.symbols('a b c h m')

dS_dt = S * (h - m * S / F)
dF_dt = F * (a - b * F - c * S)

# 2. Substitute the given parameter values
params = {a: 2, b: 1, c: 1, h: 1, m: 1}
G = dS_dt.subs(params)
H = dF_dt.subs(params)

# 3. Find the non-trivial equilibrium point (S0, F0)
# We solve the system G=0, H=0 for S>0, F>0
sols = sympy.solve([G, H], (S, F), dict=True)
non_trivial_eq = [sol for sol in sols if sol[S] > 0 and sol[F] > 0][0]
S0 = non_trivial_eq[S]
F0 = non_trivial_eq[F]

# 4. Compute the Jacobian matrix of the system (G, H) with respect to (S, F)
J = sympy.Matrix([G, H]).jacobian(sympy.Matrix([S, F]))

# 5. Evaluate the Jacobian at the equilibrium point (S0, F0)
J_eq = J.subs({S: S0, F: F0})

# The coefficients a_ij are the elements of the evaluated Jacobian
a11 = J_eq[0, 0]
a12 = J_eq[0, 1]
a21 = J_eq[1, 0]
a22 = J_eq[1, 1]

# 6. The b_ij terms are zero for a linearization at an equilibrium point
b11 = 0
b22 = 0

# 7. Print the required coefficients and the final equation
print("The coefficients for the linearized system are:")
print(f"a11 = {a11}")
print(f"a12 = {a12}")
print(f"a21 = {a21}")
print(f"a22 = {a22}")
print(f"b11 = {b11}")
print(f"b22 = {b22}")

print("\nThe final equation with these coefficients is:")
print("| x'(t) |   | {:2}  {:2} | | x(t) |   | {} |".format(int(a11), int(a12), int(b11)))
print("| y'(t) | = | {:2}  {:2} | | y(t) | + | {} |".format(int(a21), int(a22), int(b22)))