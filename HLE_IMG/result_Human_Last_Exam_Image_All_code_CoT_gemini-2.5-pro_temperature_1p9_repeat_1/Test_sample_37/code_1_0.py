import numpy as np

# --- 1. Define given values ---
# Permeability constant
mu0_over_4pi = 1e-7  # T*m/A

# Field point (target) in Cartesian coordinates
xf, yf, zf = 0.1, -0.1, 0.2
r_f = np.array([xf, yf, zf])

# Source point in Cylindrical coordinates
r_s_cyl = 0.34   # m
theta_s_deg = 60 # degrees
z_s_cyl = 0.5    # m

# Azimuthal current density magnitude
J_theta_mag = 40  # A/m^2

# --- 2. Coordinate Transformations ---
# Convert source angle to radians
theta_s_rad = np.deg2rad(theta_s_deg)

# Convert source point to Cartesian coordinates
xs = r_s_cyl * np.cos(theta_s_rad)
ys = r_s_cyl * np.sin(theta_s_rad)
zs = z_s_cyl
r_s = np.array([xs, ys, zs])

# Convert the azimuthal current density vector to Cartesian coordinates
# The direction is given by the azimuthal unit vector theta_hat = [-sin(theta), cos(theta), 0]
theta_hat = np.array([-np.sin(theta_s_rad), np.cos(theta_s_rad), 0])
J_vec = J_theta_mag * theta_hat

# --- 3. Vector Calculations ---
# Calculate the displacement vector R
R_vec = r_f - r_s
Rx, Ry, Rz = R_vec[0], R_vec[1], R_vec[2]

# Calculate the magnitude of R
R_mag = np.linalg.norm(R_vec)
R3 = R_mag**3

# Calculate the z-component of the cross product J x R
# (J x R)_z = Jx*Ry - Jy*Rx
J_cross_R_z = J_vec[0] * R_vec[1] - J_vec[1] * R_vec[0]

# --- 4. Final Calculation and Output ---
# Calculate the longitudinal component of the magnetic field, Bz
Bz = mu0_over_4pi * J_cross_R_z / R3

# Print the breakdown of the calculation as requested
print("The longitudinal component of the magnetic field (Bz) is calculated using the Biot-Savart Law:")
print("Bz = (μ₀/4π) * (J × R)z / R³\n")

print(f"Given values:")
print(f"  Field point r_f = ({xf}, {yf}, {zf}) m")
print(f"  Source point (cylindrical) = ({r_s_cyl} m, {theta_s_deg}°, {z_s_cyl} m)")
print(f"  Current density Jθ = {J_theta_mag} A/m²")
print(f"  μ₀/4π = {mu0_over_4pi} T·m/A\n")

print("Step 1: Convert source coordinates to Cartesian")
print(f"  Source point r_s = ({r_s[0]:.4f}, {r_s[1]:.4f}, {r_s[2]:.4f}) m")
print(f"  Current density vector J = ({J_vec[0]:.4f}, {J_vec[1]:.4f}, {J_vec[2]:.4f}) A/m²\n")

print("Step 2: Calculate the displacement vector R = r_f - r_s")
print(f"  R = ({Rx:.4f}, {Ry:.4f}, {Rz:.4f}) m\n")

print("Step 3: Calculate the terms for the Biot-Savart Law")
print(f"  Magnitude R = |R| = {R_mag:.4f} m")
print(f"  R³ = {R3:.4f} m³")
print(f"  (J × R)z = Jx*Ry - Jy*Rx = ({J_vec[0]:.4f})*({Ry:.4f}) - ({J_vec[1]:.4f})*({Rx:.4f}) = {J_cross_R_z:.4f} A·m/m²\n")

print("Step 4: Calculate the final result for Bz")
print(f"Bz = {mu0_over_4pi} * ({J_cross_R_z:.4f}) / {R3:.4f}")
print(f"Bz = {Bz:.4e} T")

final_answer = Bz
<<<1.2016538462725832e-05>>>