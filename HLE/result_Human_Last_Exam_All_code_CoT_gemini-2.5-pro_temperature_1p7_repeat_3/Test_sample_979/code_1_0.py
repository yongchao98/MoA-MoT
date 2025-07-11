import sympy as sp

# Define symbols
H, r, theta, R, K0, mu, mu0 = sp.symbols('H r theta R K_0 mu mu_0', real=True, positive=True)
z_hat, r_hat, theta_hat = sp.symbols('z_hat r_hat theta_hat')

# The magnetic field inside the sphere (r < R)
# Based on our derivation, this corresponds to the expression in option E
H_in_coeff = (2 * mu0 / mu) * (K0 / (1 + 2 * mu0 / mu))
H_in_expr = H_in_coeff * z_hat

# The magnetic field outside the sphere (r > R)
# Based on our derivation, this corresponds to the expression in option E
H_out_coeff = (K0 / (1 + 2 * mu0 / mu)) * (R**3 / r**3)
H_out_vector = (2 * sp.cos(theta) * r_hat + sp.sin(theta) * theta_hat)
H_out_expr = H_out_coeff * H_out_vector

# Print the final result corresponding to choice E
print("The correct answer is E.")
print("The magnetic field H(r, theta) is given by:")
print("\nFor 0 < r < R:")
print(f"  H_in = ( (2 * mu_0) / mu ) * ( K_0 / (1 + (2 * mu_0)/mu) ) * z_hat")
# A simplified but equivalent form from derivation
print(f"       = ( (2 * mu_0 * K_0) / (mu + 2 * mu_0) ) * z_hat")

print("\nFor R < r < infinity:")
print(f"  H_out = ( K_0 / (1 + (2 * mu_0)/mu) ) * (R^3 / r^3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)")
# A simplified but equivalent form from derivation
print(f"        = ( (mu * K_0 * R^3) / ((mu + 2*mu_0) * r^3) ) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)")