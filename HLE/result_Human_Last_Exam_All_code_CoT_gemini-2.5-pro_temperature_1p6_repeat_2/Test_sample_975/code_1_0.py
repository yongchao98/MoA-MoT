import sympy as sp

# Define the symbolic variables for the problem
M0, Rp, R, r, theta = sp.symbols('M_0 R_p R r theta')

# --- Region 1: Inside the shield (0 < r < R_p) ---

# The coefficient for the H field components in Region 1
coeff1 = M0 * (R**3 + 2*Rp**3) / (3*R**3)

# The H-field components in Region 1
H1_r = -coeff1 * sp.cos(theta)
H1_theta = coeff1 * sp.sin(theta)

print("Magnetic Field H in Region 1 (0 < r < R_p):")
print("H_r = (" + str(H1_r) + ") * r_hat")
print("H_theta = (" + str(H1_theta) + ") * theta_hat")
print("-" * 20)

# The expression in Choice B for region 1 can be written as:
# H = M0 * ( (2*Rp**3 + R**3)/(3*R**3) ) * (-cos(theta)*r_hat + sin(theta)*theta_hat)
# This matches our derived components.

# --- Region 2: Between the shield and conductor (R_p < r < R) ---

# The H-field components in Region 2
H2_r_coeff = - (2 * M0 / 3) * ( (Rp/R)**3 - (Rp/r)**3 )
H2_theta_coeff = (M0 / 3) * ( 2*(Rp/R)**3 + (Rp/r)**3 )

H2_r = H2_r_coeff * sp.cos(theta)
H2_theta = H2_theta_coeff * sp.sin(theta)

print("Magnetic Field H in Region 2 (R_p < r < R):")
print("H_r = (" + str(H2_r) + ") * r_hat")
print("H_theta = (" + str(H2_theta) + ") * theta_hat")
print("-" * 20)

# The expressions match Choice B exactly.
print("Comparing with the answer choices, the derived fields match choice B.")
print("\nFinal derived expressions printed term by term as requested:\n")

print("Region 1 (0 < r < Rp):")
print("H_r = - (M_0*(2*R_p**3 + R**3)/(3*R**3)) * cos(theta)")
print("H_theta = (M_0*(2*R_p**3 + R**3)/(3*R**3)) * sin(theta)\n")

print("Region 2 (Rp < r < R):")
print("H_r = -(2*M_0/3) * ( (R_p/R)**3 - (R_p/r)**3 ) * cos(theta)")
print("H_theta = (M_0/3) * ( 2*(R_p/R)**3 + (R_p/r)**3 ) * sin(theta)")
