import sympy as sp

# Define symbolic variables
# Note: mu_0 is represented as mu0 in the code for simplicity.
K0, R, r, theta, mu, mu0 = sp.symbols('K_0 R r theta mu mu_0')

# Common denominator term
denominator = (1 + 2 * mu0 / mu)

# Field inside the sphere (0 < r < R)
# This is a uniform field in the z-direction.
H_in_coeff = (K0 * (2 * mu0 / mu)) / denominator
H_in_direction = "z_hat (unit vector in z-direction)"

# Field outside the sphere (R < r < infinity)
# This is a dipole field.
H_out_coeff = (K0 / denominator) * (R**3 / r**3)
H_out_radial_comp = 2 * sp.cos(theta)
H_out_polar_comp = sp.sin(theta)

print("The magnetic field H(r, theta) is determined in two regions:")
print("-" * 50)

# Print the solution for the region inside the sphere
print("For the region inside the sphere (0 < r < R):")
print("H_in = (C_in) * z_hat")
print(f"The coefficient C_in is:")
sp.pprint(H_in_coeff, use_unicode=True)
print("So, the full expression is:")
print(f"H_in(r, theta) = ({sp.pretty(H_in_coeff)}) * {H_in_direction}\n")


# Print the solution for the region outside the sphere
print("For the region outside the sphere (R < r < infinity):")
print("H_out = (C_out) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)")
print(f"The coefficient C_out is:")
sp.pprint(H_out_coeff, use_unicode=True)
print("So, the full expression is:")
print(f"H_out(r, theta) = ({sp.pretty(H_out_coeff)}) * ({sp.pretty(H_out_radial_comp)}*r_hat + {sp.pretty(H_out_polar_comp)}*theta_hat)\n")

print("-" * 50)
print("This result corresponds to option E.")