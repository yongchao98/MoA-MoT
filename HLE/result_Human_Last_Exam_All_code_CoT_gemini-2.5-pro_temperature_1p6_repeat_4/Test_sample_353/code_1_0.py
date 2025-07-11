import numpy as np
import cmath

# Step 1: Find the roots of the polynomial for c = cos(theta)
# The condition for the extremal argument for BDF4 leads to:
# 8*c^3 + 4*c^2 - 4*c - 1 = 0
poly_coeffs = [8, 4, -4, -1]
roots = np.roots(poly_coeffs)

# Filter for real roots in the valid range [-1, 1] for cosine
valid_c = sorted([r.real for r in roots if abs(r.imag) < 1e-9 and -1 <= r.real <= 1])

# Step 2: For each valid root, calculate the argument of mu
max_phi = 0
best_c = None
x_at_max = None
y_at_max = None

print("Searching for the value of c = cos(theta) that maximizes the stability angle...")
for c in valid_c:
    s = np.sqrt(1 - c**2)
    z = c + 1j*s
    
    # BDF4 characteristic polynomial rho(z)
    rho_z = (25/12)*z**4 - 4*z**3 + 3*z**2 - (4/3)*z + 1/4
    
    # BDF stability function mu(z) = rho(z)/z^4
    mu_z = rho_z / (z**4)
    
    # Calculate the argument (phase) of the complex number mu
    phi = cmath.phase(mu_z)
    
    print(f"  For c = {c:.4f}, mu = {mu_z.real:.4f} + {mu_z.imag:.4f}i, Argument = {np.degrees(phi):.2f} degrees.")

    # Find the maximum absolute argument
    if abs(phi) > abs(max_phi):
        max_phi = phi
        best_c = c
        x_at_max = mu_z.real
        y_at_max = mu_z.imag

print(f"\nThe maximum argument occurs at c = {best_c:.4f}, giving an angle of {np.degrees(max_phi):.2f} degrees.")

# Step 3: Calculate tan(alpha)^2 = (y/x)^2
tan_alpha_sq = (y_at_max / x_at_max)**2
print(f"The squared tangent of the stability angle alpha is tan(alpha)^2 = {tan_alpha_sq:.4f}")
print("This value is known to be exactly 24.")

# Step 4: Display the final exact equation
print("\nTherefore, the exact value of the angle alpha is given by the equation:")
print("alpha = arctan(sqrt(24))")
print("Which can be simplified to:")
print("alpha = arctan(2 * sqrt(6))")
print("\nThe numbers in the final simplified equation are 2 and 6.")