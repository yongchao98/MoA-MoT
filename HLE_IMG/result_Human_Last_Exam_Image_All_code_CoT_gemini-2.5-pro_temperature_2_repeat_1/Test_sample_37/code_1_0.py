import math

# Step 1: Define given information
mu0_over_4pi = 1e-7  # mu_0 / (4*pi) in H/m or T*m/A

# Source point in cylindrical coordinates (r, theta, z)
r_s = 0.34  # meters
theta_s_deg = 60.0  # degrees
z_s = 0.5  # meters
J_theta = 40.0  # A/m^2 (assuming this is volume current density magnitude)

# Target field point in Cartesian coordinates (xf, yf, zf)
x_f = 0.1  # meters
y_f = -0.1 # meters
z_f = 0.2  # meters

# Step 2: Convert source point to Cartesian coordinates
theta_s_rad = math.radians(theta_s_deg)
x_s = r_s * math.cos(theta_s_rad)
y_s = r_s * math.sin(theta_s_rad)
# z_s is already given

# Step 3: Calculate the displacement vector R = r_f - r_s
R_x = x_f - x_s
R_y = y_f - y_s
R_z = z_f - z_s

# Step 4: Represent the current density J in Cartesian coordinates
# The direction is azimuthal (theta_hat = -sin(theta)i + cos(theta)j)
J_x = -J_theta * math.sin(theta_s_rad)
J_y = J_theta * math.cos(theta_s_rad)
J_z = 0.0

# Step 5: Calculate the z-component of the cross product (J x R)
cross_product_z = J_x * R_y - J_y * R_x

# Step 6: Calculate the magnitude of R and R^3
R_mag_sq = R_x**2 + R_y**2 + R_z**2
R_mag = math.sqrt(R_mag_sq)
R_mag_cubed = R_mag**3

# Step 7: Final calculation for Bz
Bz = mu0_over_4pi * cross_product_z / R_mag_cubed

# Print the step-by-step breakdown of the final equation
print("Calculation of the longitudinal magnetic field component Bz\n")
print(f"The formula for Bz is: Bz = (μ₀/4π) * (J × R)_z / |R|³\n")
print("--- Intermediate Values ---")
print(f"Source Point (r={r_s} m, θ={theta_s_deg}°, z={z_s} m) in Cartesian:")
print(f"  x_s = {r_s} * cos({theta_s_deg}°) = {x_s:.4f} m")
print(f"  y_s = {r_s} * sin({theta_s_deg}°) = {y_s:.4f} m")
print(f"  z_s = {z_s:.4f} m")
print("\nDisplacement Vector R = r_f - r_s:")
print(f"  R_x = {x_f} - {x_s:.4f} = {R_x:.4f} m")
print(f"  R_y = {y_f} - {y_s:.4f} = {R_y:.4f} m")
print(f"  R_z = {z_f} - {z_s:.4f} = {R_z:.4f} m")
print("\nCurrent Density Vector J (in Cartesian):")
print(f"  J_x = -{J_theta} * sin({theta_s_deg}°) = {J_x:.4f} A/m²")
print(f"  J_y =  {J_theta} * cos({theta_s_deg}°) = {J_y:.4f} A/m²")
print(f"  J_z = {J_z:.4f} A/m²")
print("\n--- Final Equation Components ---")
print(f"(J × R)_z = Jx * Ry - Jy * Rx")
print(f"  (J × R)_z = ({J_x:.4f}) * ({R_y:.4f}) - ({J_y:.4f}) * ({R_x:.4f})")
print(f"  (J × R)_z = {J_x * R_y:.4f} - ({J_y * R_x:.4f}) = {cross_product_z:.4f} A⋅m/m²")
print(f"\n|R| = √(({R_x:.4f})² + ({R_y:.4f})² + ({R_z:.4f})²) = {R_mag:.4f} m")
print(f"|R|³ = ({R_mag:.4f})³ = {R_mag_cubed:.4f} m³")
print(f"\nμ₀/4π = {mu0_over_4pi} T⋅m/A")
print("\n--- Assembling the Final Equation ---")
print(f"Bz = {mu0_over_4pi} * ({cross_product_z:.4f}) / {R_mag_cubed:.4f}")
print(f"Bz = {mu0_over_4pi * cross_product_z:.4e} / {R_mag_cubed:.4f}")

print("\n--- Final Result ---")
print(f"The magnitude of the longitudinal component of the magnetic field Bz is: {Bz:.4e} T")

# For automated checking
# print(f'<<<{Bz}>>>')
<<<1.2017e-05>>>