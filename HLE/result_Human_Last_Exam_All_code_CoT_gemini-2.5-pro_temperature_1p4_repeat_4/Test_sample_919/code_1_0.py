import sympy as sp

# Define symbols for the variables in the problem
mu_0, mu, K_0, a, y, d = sp.symbols('mu_0 mu K_0 a y d', real=True, positive=True)

# The derived magnetic field intensity H at the conductor surface (x=d)
# H_y is the only non-zero component at x=d.
Denominator_H = sp.cosh(a*d) + (mu_0/mu) * sp.sinh(a*d)
H_y_at_d = (K_0 * sp.sin(a*y)) / Denominator_H

# The magnitude squared of the magnetic field at x=d
H_magnitude_squared_at_d = H_y_at_d**2

# The formula for the force per unit area on the conductor
# f/A = (mu_0 / 2) * |H|^2, in the x-direction
force_per_area_x = (mu_0 / 2) * H_magnitude_squared_at_d

# Print the final formula for the force per unit area in a readable format
# The formula corresponds to option C.
print("The force per unit area on the x=d interface is:")
# Manually format the string for clarity, as sympy's pretty print can be complex
numerator_str = f"mu_0 * K_0**2 * sin(a*y)**2"
denominator_str = f"2 * [cosh(a*d) + (mu_0/mu) * sinh(a*d)]**2"
print(f"f/area = ({numerator_str}) / ({denominator_str}) * i_x")

print("\nComparing with the options, this matches option C:")
print(f"C. f/area = (mu_0/2) * (K_0**2 * sin(a*y)**2) / ([cosh(a*d) + (mu_0/mu) * sinh(a*d)]**2) * i_x")