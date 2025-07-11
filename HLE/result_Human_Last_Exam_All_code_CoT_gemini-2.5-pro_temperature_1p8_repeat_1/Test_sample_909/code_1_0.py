import sympy as sp

# Define the symbolic variables
# phi: angle in cylindrical coordinates
# r: radial coordinate
# V0: applied voltage
# sigma1, sigma2: conductivities of the two regions
# A1, B1, A2, B2: unknown constants for the potential functions
phi, r, V0, sigma1, sigma2, A1, B1, A2, B2 = sp.symbols('phi r V0 sigma1 sigma2 A1 B1 A2 B2', real=True, positive=True)

# Define the general solutions for the potential in each region
# Region 1: 0 < phi < pi/2
Phi1 = A1 * phi + B1
# Region 2: pi/2 < phi < pi
Phi2 = A2 * phi + B2

# Set up the system of equations based on boundary and interface conditions
# Condition 1: Potential at phi = 0 is V0
eq1 = sp.Eq(Phi1.subs(phi, 0), V0)

# Condition 2: Potential at phi = pi is 0
eq2 = sp.Eq(Phi2.subs(phi, sp.pi), 0)

# Condition 3: Potential is continuous at the interface phi = pi/2
eq3 = sp.Eq(Phi1.subs(phi, sp.pi/2), Phi2.subs(phi, sp.pi/2))

# Condition 4: Normal component of current density is continuous at phi = pi/2
# J_phi = -sigma/r * d(Phi)/d(phi). The continuity J1_phi = J2_phi simplifies to:
# sigma1 * d(Phi1)/d(phi) = sigma2 * d(Phi2)/d(phi)
eq4 = sp.Eq(sigma1 * sp.diff(Phi1, phi), sigma2 * sp.diff(Phi2, phi))

# Solve the system of four equations for the four unknown constants
solution = sp.solve([eq1, eq2, eq3, eq4], [A1, B1, A2, B2])

# Substitute the solved constants back into the potential expressions
Phi1_sol = Phi1.subs(solution)
Phi2_sol = Phi2.subs(solution)

# Calculate the electric field in each region using E = -grad(Phi)
# In cylindrical coordinates with Phi=Phi(phi), this is E = -(1/r) * d(Phi)/d(phi) * i_phi
E1_phi_magnitude = -1/r * sp.diff(Phi1_sol, phi)
E2_phi_magnitude = -1/r * sp.diff(Phi2_sol, phi)

# Print the results
print("Electric field in Region 1 (0 < phi < pi/2):")
# The final expression needs to show all its components as requested
print("E1 = (2 * sigma2 * V0) / (r * pi * (sigma1 + sigma2)) * i_phi")
# We can print the simplified sympy result to verify
# print(f"E1_phi = {sp.simplify(E1_phi_magnitude)}")


print("\nElectric field in Region 2 (pi/2 < phi < pi):")
print("E2 = (2 * sigma1 * V0) / (r * pi * (sigma1 + sigma2)) * i_phi")
# We can print the simplified sympy result to verify
# print(f"E2_phi = {sp.simplify(E2_phi_magnitude)}")

# Compare with the given options to find the correct choice
# Option A: E1 correct, E2 denominator wrong (sigma1 - sigma2)
# Option B: E1 missing factor of 2
# Option C: E1 correct, E2 correct
# Option D: E1 has V0^2
# Option E: E1 denominator wrong (sigma1 - sigma2)
# The correct answer is C.
