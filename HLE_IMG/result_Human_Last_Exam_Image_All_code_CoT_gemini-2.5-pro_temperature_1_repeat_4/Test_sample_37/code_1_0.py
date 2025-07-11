import math

# --- Given Information & Constants ---
# Permeability of free space constant
mu0_div_4pi = 1e-7  # H/m

# Target field point (Cartesian)
xf, yf, zf = 0.1, -0.1, 0.2  # m

# Source point (Cylindrical)
r_s = 0.34  # m
theta_s_deg = 60  # degrees
z_s = 0.5  # m

# Azimuthal surface current density (assuming A/m)
K_theta = 40.0  # A/m

# --- Step-by-Step Calculation ---

# Step 1: Convert source point to Cartesian coordinates (xs, ys, zs)
# Convert theta from degrees to radians for math functions
theta_s_rad = math.radians(theta_s_deg)
xs = r_s * math.cos(theta_s_rad)
ys = r_s * math.sin(theta_s_rad)
# zs is the same in both coordinate systems

print(f"Step 1: Convert source point to Cartesian coordinates")
print(f"Source point (r, θ, z) = ({r_s} m, {theta_s_deg}°, {z_s} m)")
print(f"Source point (x_s, y_s, z_s) = ({xs:.4f} m, {ys:.4f} m, {z_s:.4f} m)\n")

# Step 2: Calculate the displacement vector R = rf - rs
Rx = xf - xs
Ry = yf - ys
Rz = zf - z_s
R_mag = math.sqrt(Rx**2 + Ry**2 + Rz**2)

print(f"Step 2: Calculate the displacement vector R = r_f - r_s")
print(f"Target point (x_f, y_f, z_f) = ({xf} m, {yf} m, {zf} m)")
print(f"Displacement vector R = ({Rx:.4f}, {Ry:.4f}, {Rz:.4f}) m")
print(f"Magnitude |R| = {R_mag:.4f} m\n")

# Step 3: Represent the surface current density vector K in Cartesian coordinates
# K is in the azimuthal direction (theta-hat), which is -sin(θ)i + cos(θ)j
Kx = -K_theta * math.sin(theta_s_rad)
Ky = K_theta * math.cos(theta_s_rad)

print(f"Step 3: Represent the current density vector K in Cartesian coordinates")
print(f"Surface current density K_theta = {K_theta} A/m")
print(f"Current vector K = ({Kx:.4f}, {Ky:.4f}, 0) A/m\n")

# Step 4: Calculate the z-component of the cross product (K x R)
K_cross_R_z = Kx * Ry - Ky * Rx

print(f"Step 4: Calculate the z-component of the cross product (K x R)")
print(f"(K x R)_z = Kx * Ry - Ky * Rx")
print(f"(K x R)_z = ({Kx:.4f}) * ({Ry:.4f}) - ({Ky:.4f}) * ({Rx:.4f})")
print(f"(K x R)_z = {K_cross_R_z:.4f} A*m\n")

# Step 5: Calculate Bz using the Biot-Savart Law
# Bz = (μ₀/4π) * (K × R)_z / |R|³
Bz = mu0_div_4pi * K_cross_R_z / (R_mag**3)

print(f"Step 5: Calculate the final magnetic field component B_z")
print(f"B_z = (μ₀/4π) * (K × R)_z / |R|³")
print(f"B_z = ({mu0_div_4pi}) * ({K_cross_R_z:.4f}) / ({R_mag:.4f})³")
print(f"B_z = ({mu0_div_4pi}) * ({K_cross_R_z:.4f}) / ({R_mag**3:.4f})")
print(f"The magnitude of the longitudinal component of the magnetic field is:")
print(f"B_z = {Bz:.4e} T")
