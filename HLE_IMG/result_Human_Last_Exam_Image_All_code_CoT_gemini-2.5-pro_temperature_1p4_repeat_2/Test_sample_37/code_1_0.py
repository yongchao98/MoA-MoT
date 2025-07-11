import math

# Step 1: Define constants and given values
mu0_over_4pi = 1e-7  # Permeability of free space constant in H/m

# Target field point (Cartesian)
xf, yf, zf = 0.1, -0.1, 0.2

# Source point (Cylindrical)
r_s_cyl = 0.34
theta_s_deg = 60
z_s_cyl = 0.5

# Azimuthal current density (Assuming A/m)
K_theta = 40.0

print("--- Problem Setup ---")
print(f"Target field point (xf, yf, zf): ({xf}, {yf}, {zf}) m")
print(f"Source point (r, θ, z): ({r_s_cyl} m, {theta_s_deg}°, {z_s_cyl} m)")
print(f"Azimuthal current density K_theta: {K_theta} A/m\n")

# Step 2: Convert source point to Cartesian coordinates
theta_s_rad = math.radians(theta_s_deg)
xs = r_s_cyl * math.cos(theta_s_rad)
ys = r_s_cyl * math.sin(theta_s_rad)
zs = z_s_cyl
print("--- Step 1: Coordinate Conversion ---")
print(f"Source point in Cartesian (xs, ys, zs): ({xs:.5f}, {ys:.5f}, {zs:.5f}) m")


# Step 3: Calculate the displacement vector R = r_f - r_s
Rx = xf - xs
Ry = yf - ys
Rz = zf - zs
print(f"Displacement vector R = (Rx, Ry, Rz): ({Rx:.5f}, {Ry:.5f}, {Rz:.5f}) m\n")

# Step 4: Convert current density vector K to Cartesian coordinates
# K = K_theta * theta_hat, where theta_hat = -sin(theta)x_hat + cos(theta)y_hat
Kx = K_theta * (-math.sin(theta_s_rad))
Ky = K_theta * (math.cos(theta_s_rad))
Kz = 0.0 # No z-component for azimuthal current
print("--- Step 2: Vector Calculation ---")
print(f"Current density vector K in Cartesian (Kx, Ky, Kz): ({Kx:.5f}, {Ky:.5f}, {Kz:.1f}) A/m")

# Step 5: Calculate the z-component of the cross product (K x R)_z
cross_product_z = Kx * Ry - Ky * Rx
print(f"Z-component of cross product (K x R)z: {cross_product_z:.5f} A")

# Step 6: Calculate the magnitude of R cubed (R^3)
R_sq = Rx**2 + Ry**2 + Rz**2
R = math.sqrt(R_sq)
R_cubed = R**3
print(f"Magnitude of R: {R:.5f} m")
print(f"Magnitude of R-cubed: {R_cubed:.5f} m^3\n")

# Step 7: Calculate the longitudinal component of the magnetic field (Bz)
Bz = mu0_over_4pi * cross_product_z / R_cubed

# Final Output
print("--- Final Calculation ---")
print("The equation for the z-component of the magnetic field is:")
print("Bz = (μ₀/4π) * (Kx * Ry - Ky * Rx) / (Rx² + Ry² + Rz²)^(3/2)\n")

print("Substituting the calculated values into the equation:")
print(f"Bz = {mu0_over_4pi} * (({Kx:.5f}) * ({Ry:.5f}) - ({Ky:.5f}) * ({Rx:.5f})) / (({Rx:.5f})² + ({Ry:.5f})² + ({Rz:.5f})²)^(3/2)")
print(f"Bz = {mu0_over_4pi} * ({cross_product_z:.5f}) / ({R_cubed:.5f})")
print(f"Bz = {mu0_over_4pi} * {cross_product_z/R_cubed:.5f}")

print("\nThe final magnitude of the longitudinal component of the magnetic field is:")
print(f"Bz = {Bz:.5e} T")
