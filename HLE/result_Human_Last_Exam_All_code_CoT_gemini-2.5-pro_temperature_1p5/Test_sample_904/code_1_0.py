# Goal: Programmatically display the derived coefficients and the governing equation.

# The derived coefficients for the linear governing equation
# A(r) * d^2(xi)/dr^2 + B(r) * d(xi)/dr + C(r, xi) = 0

# From the derivation, A(r) is the coefficient for the second derivative term.
# It arises directly from the surface tension 'gamma'.
A_r = "gamma"

# B(r) is the coefficient for the first derivative term.
# It comes from the surface tension 'gamma' and the radial 'r' term in the
# linearized curvature expression for cylindrical coordinates.
B_r = "gamma / r"

# Print the explanation and results
print("The governing linear equation for the interfacial shape xi(r) is derived by balancing the capillary pressure with external forces.")
print("The general form of the equation is:")
print("  A(r) * d^2(xi)/dr^2 + B(r) * d(xi)/dr + C(r, xi) = 0\n")

print("From the derivation based on the linearized Young-Laplace equation in cylindrical coordinates, the coefficients A(r) and B(r) are identified as:")
print(f"  A(r) = {A_r}")
print(f"  B(r) = {B_r}")

print("\nWhere:")
print("  'gamma' is the surface tension between the two fluids.")
print("  'r' is the radial coordinate.")
print("  'C(r, xi)' is a term representing all other forces, primarily the electrostatic pressure in this system.")

print("\nThe final equation with each term shown is:")
print(f"({A_r}) * d^2(xi)/dr^2 + ({B_r}) * d(xi)/dr + C(r, xi) = 0")