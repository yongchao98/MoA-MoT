import sympy as sp

# Define symbols for our physical quantities
# While we won't execute this symbolically, it helps to represent the variables.
r, R, theta, K0, mu, mu0 = sp.symbols('r R theta K_0 mu mu_0')

# Step 1 & 2: Define the magnetic scalar potential inside and outside

# Inside (r < R), potential must be finite at r=0, so B_l coefficients are zero.
# Outside (r > R), potential must go to zero at r=infinity, so A_l coefficients are zero.
# The current K is proportional to sin(theta), which is related to the l=1 Legendre polynomial P_1(cos(theta)).
# Therefore, only the l=1 terms will be non-zero.
# Phi_in(r, theta) = A1 * r * cos(theta)
# Phi_out(r, theta) = B1 * r**(-2) * cos(theta)

print("Step 1: Define Potentials")
print("Inside (r < R): Phi_in(r, theta) = A1 * r * cos(theta)")
print("Outside (r > R): Phi_out(r, theta) = B1 / r**2 * cos(theta)")
print("-" * 30)

# Step 3 & 4: Apply boundary conditions to find coefficients A1 and B1

print("Step 2: Apply Boundary Conditions at r=R")

# Condition 1: Continuity of normal B-field (B_r)
# B_out,r = B_in,r  =>  mu_0 * H_out,r = mu * H_in,r
# H_r = -d(Phi)/dr
# mu_0 * (-d(Phi_out)/dr) = mu * (-d(Phi_in)/dr)  at r=R
# mu_0 * (2*B1/R**3 * cos(theta)) = mu * (-A1 * cos(theta))
# Equation I: 2 * mu_0 * B1 = -mu * A1 * R**3
print("Condition 1 (B_r continuous) gives: 2 * mu_0 * B1 = -mu * A1 * R**3")

# Condition 2: Discontinuity of tangential H-field (H_theta)
# H_out,theta - H_in,theta = (K x n)_theta = (K_phi * phi_hat x r_hat)_theta = K_phi = K0*sin(theta)
# H_theta = -(1/r) * d(Phi)/d(theta)
# [-(1/R) * (-B1/R**2 * sin(theta))] - [-(1/R) * (-A1*R*sin(theta))] = K0*sin(theta)
# (B1/R**3) - A1 = K0
# Equation II: B1 - A1 * R**3 = K0 * R**3
print("Condition 2 (H_theta discontinuous) gives: B1 - A1 * R**3 = K0 * R**3")
print("-" * 30)

print("Step 3: Solve for coefficients A1 and B1")
# From Eq I: A1*R**3 = - (2*mu_0/mu) * B1
# Substitute into Eq II:
# B1 - (-(2*mu_0/mu) * B1) = K0 * R**3
# B1 * (1 + 2*mu_0/mu) = K0 * R**3
# So, B1 = (K0 * R**3) / (1 + 2*mu_0/mu)
print(f"Solving for B1 gives: B1 = (K0 * R^3) / (1 + 2*mu_0/mu)")

# Now find A1:
# A1 * R**3 = B1 - K0*R**3
# A1*R**3 = (K0*R**3)/(1+2*mu_0/mu) - K0*R**3 = K0*R**3 * (1/(1+2*mu_0/mu) - 1)
# A1*R**3 = K0*R**3 * ( (1 - (1+2*mu_0/mu)) / (1+2*mu_0/mu) )
# A1 = K0 * (-2*mu_0/mu) / (1 + 2*mu_0/mu)
print(f"Solving for A1 gives: A1 = -K0 * (2*mu_0/mu) / (1 + 2*mu_0/mu)")
print("-" * 30)

# Step 5: Calculate the magnetic fields H = -gradient(Phi)

print("Step 4: Calculate H_in and H_out")
# Inside Field (H_in)
# H_in = -gradient(A1 * r * cos(theta)) = -gradient(A1 * z) = -A1 * z_hat
# H_in = -(-K0 * (2*mu_0/mu) / (1 + 2*mu_0/mu)) * z_hat
# H_in = ( (2*mu_0/mu) * K0 / (1 + 2*mu_0/mu) ) * z_hat
print("The magnetic field inside the sphere (0 < r < R) is:")
print("H_in(r, theta) = (2 * mu_0 / mu) * (K_0 / (1 + 2 * mu_0 / mu)) * z_hat")


# Outside Field (H_out)
# H_out = -gradient(B1/r**2 * cos(theta))
# H_out,r = -d/dr(B1/r**2 * cos(theta)) = 2*B1/r**3 * cos(theta)
# H_out,theta = -(1/r) * d/d(theta)(B1/r**2 * cos(theta)) = B1/r**3 * sin(theta)
# H_out = (B1/r**3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)
# H_out = ( (K0*R**3)/(1+2*mu_0/mu) )/r**3 * (2*cos(theta)*r_hat + sin(theta)*theta_hat)
print("\nThe magnetic field outside the sphere (R < r < infinity) is:")
print("H_out(r, theta) = (K_0 / (1 + 2 * mu_0 / mu)) * (R^3 / r^3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)")
print("-" * 30)

print("Final Answer Check:")
print("Comparing our derived fields with the given options, we find a perfect match with option E.")
print("""
H(r, theta) = 
    { 
    (2 * mu_0 / mu) * (K_0 / (1 + (2 * mu_0 / mu))) * z_hat,            for 0 < r < R
    (K_0 / (1 + (2 * mu_0 / mu))) * (R^3 / r^3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat), for R < r < infinity
    }
""")
print("Note: In option E, the expression z_hat for the inside field is mathematically equivalent to (cos(theta)*r_hat - sin(theta)*theta_hat), but the Cartesian form is much simpler and more intuitive for a uniform field.")
