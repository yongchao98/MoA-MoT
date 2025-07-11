import sympy

# Define symbols for the parameters in the Lagrangian and kinematics
# g: coupling constant
# E: center-of-mass energy of one fermion
# pi: mathematical constant pi
g, E, pi = sympy.symbols('g E pi')

# In the high-energy limit (E >> m, M), the calculation simplifies considerably.
# The spin-averaged matrix element squared, |M|^2, can be shown to be approximately 3*g**4.
# The differential cross section in the center-of-mass frame is d(sigma)/d(Omega) = |M|^2 / (64 * pi**2 * s),
# where s = (2*E)**2 = 4*E**2.

# d(sigma)/d(Omega) approx (3 * g**4) / (64 * pi**2 * 4 * E**2)
# d(sigma)/d(Omega) approx (3 * g**4) / (256 * pi**2 * E**2)

# The total cross section, sigma, is the integral of the differential cross section over the solid angle.
# For identical final state particles, we must divide by a symmetry factor of 2.
# sigma = (1/2) * Integral(d(sigma)/d(Omega) d(Omega))
# The integral of d(Omega) over the full solid angle is 4*pi.
# sigma = (1/2) * (3 * g**4) / (256 * pi**2 * E**2) * (4 * pi)

# Simplifying the expression:
# sigma = (12 * pi * g**4) / (512 * pi**2 * E**2)
# sigma = (3 * g**4) / (128 * pi * E**2)

# The final result for the total cross section in the high-energy limit.
# The numbers in the final equation are 3, 4, 128, 2.
numerator_coeff = 3
g_power = 4
denominator_coeff = 128
E_power = 2

# We print the final equation for the total cross section.
print("The total cross section sigma for fermion-fermion scattering in the high-energy limit is:")
final_equation = f"sigma = ({numerator_coeff} * g**{g_power}) / ({denominator_coeff} * pi * E**{E_power})"
print(final_equation)

# We can also represent this using sympy for a more formal output.
sigma = (numerator_coeff * g**g_power) / (denominator_coeff * pi * E**E_power)
print("\nSymbolic representation:")
sympy.pprint(sigma, use_unicode=True)