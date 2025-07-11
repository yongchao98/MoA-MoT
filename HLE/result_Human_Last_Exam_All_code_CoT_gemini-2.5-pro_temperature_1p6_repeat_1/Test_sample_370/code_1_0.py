import sympy

# Define all symbols used in the calculation.
# g: coupling constant
# E: center-of-mass energy of one fermion
# pi: the mathematical constant pi
g, E, pi = sympy.symbols('g E pi', real=True, positive=True)
s = sympy.Symbol('s', real=True, positive=True)

# Step 1: Explain the setup and assumptions
print("This script calculates the total cross-section for fermion-fermion scattering")
print("psi(p1) + psi(p2) -> psi(p3) + psi(p4) based on the provided Lagrangian.")
print("The calculation is performed in the high-energy limit, where E >> m, M.\n")

# Step 2: Define kinematic variables
# In the center-of-mass frame and high-energy limit, the Mandelstam variable s is:
s_val = 4 * E**2
print("The Mandelstam variable s in the CM frame is:")
print(f"s = 4 * E^2\n")

# Step 3: State the spin-averaged squared matrix element
# The calculation involves t- and u-channel diagrams and trace algebra.
# In the high-energy limit, the result simplifies to a constant.
M_bar_sq = 3 * g**4
print("The spin-averaged squared matrix element |M_bar|^2 in the high-energy limit is:")
print(f"|M_bar|^2 = 3 * g^4\n")

# Step 4: Calculate the differential cross-section
print("The differential cross-section d(sigma)/d(Omega) in the CM frame is:")
print("d(sigma)/d(Omega) = |M_bar|^2 / (64 * pi^2 * s)")
# Substitute the expressions for |M_bar|^2 and s
diff_sigma_d_omega = M_bar_sq / (64 * pi**2 * s_val)
print(f"d(sigma)/d(Omega) = (3 * g^4) / (64 * pi^2 * (4 * E^2)) = {sympy.simplify(diff_sigma_d_omega)}\n")

# Step 5: Calculate the total cross-section
print("To find the total cross-section sigma, we integrate over the solid angle (4*pi)")
print("and include a symmetry factor of 1/2 for identical final state particles.")
print("sigma = (1/2) * Integral( d(sigma)/d(Omega) d(Omega) )")

# The integrand is constant, so the integral is the integrand * 4*pi.
sigma = (sympy.Rational(1, 2)) * diff_sigma_d_omega * (4 * pi)
final_sigma = sympy.simplify(sigma)

print("\nThe step-by-step integration is:")
print(f"sigma = (1/2) * ({sympy.simplify(diff_sigma_d_omega)}) * Integral(d(Omega))")
print(f"sigma = (1/2) * ({sympy.simplify(diff_sigma_d_omega)}) * (4 * pi)")
print(f"sigma = {final_sigma}\n")


# Final result with each number explicitly shown as requested
num_coeff = 3
den_coeff = 128
print("The final equation for the total cross section is:")
print(f"sigma = ({num_coeff} * g^4) / ({den_coeff} * pi * E^2)")

# Format the final answer for easy parsing
final_expression_str = str(final_sigma)
print(f"\n<<<{final_expression_str}>>>")