import numpy as np

# --- 1. Define constants and given values ---
# Permeability of free space constant (mu_0 / 4pi)
mu_0_div_4pi = 1e-7  # Units: T*m/A or H/m

# Field point (target) in Cartesian coordinates (m)
xf, yf, zf = 0.1, -0.1, 0.2
r_f = np.array([xf, yf, zf])

# Source point in Cylindrical coordinates (m, degrees, m)
r_s_cyl, theta_s_deg, z_s = 0.34, 60.0, 0.5

# Azimuthal current density (A/m^2).
# As discussed in the plan, we interpret this as a surface current density K with magnitude 40 A/m.
# We also assume a unit surface area for the source element to calculate a finite field value.
K_theta = 40.0  # A/m
delta_S = 1.0   # m^2

# --- 2. Coordinate and Vector Calculations ---
# Convert source point angle to radians
theta_s_rad = np.deg2rad(theta_s_deg)

# Convert source point from cylindrical to Cartesian coordinates
xs = r_s_cyl * np.cos(theta_s_rad)
ys = r_s_cyl * np.sin(theta_s_rad)
# z_s is already given
r_s = np.array([xs, ys, z_s])

# Calculate the displacement vector R = r_f - r_s
R_vec = r_f - r_s
Rx, Ry, Rz = R_vec

# Calculate the magnitude of the displacement vector R
R_mag = np.linalg.norm(R_vec)

# Convert the azimuthal current density vector K to Cartesian coordinates
# The azimuthal unit vector theta_hat is (-sin(theta), cos(theta), 0)
Kx = -K_theta * np.sin(theta_s_rad)
Ky = K_theta * np.cos(theta_s_rad)
Kz = 0.0
K_vec = np.array([Kx, Ky, Kz])

# Calculate the cross product C = K x R
C_vec = np.cross(K_vec, R_vec)
Cx, Cy, Cz = C_vec

# --- 3. Calculate the longitudinal magnetic field component (Bz) ---
# Using the Biot-Savart Law: B = (mu_0 / 4pi) * (K x R) / R^3 * delta_S
Bz = mu_0_div_4pi * Cz / (R_mag**3) * delta_S

# --- 4. Print the results step-by-step ---
print("--- Calculation Steps ---")
print(f"Field Point (xf, yf, zf) = ({xf}, {yf}, {zf}) m")
print(f"Source Point (r, θ, z) = ({r_s_cyl}, {theta_s_deg}°, {z_s}) m")
print(f"Source Point in Cartesian (xs, ys, zs) = ({xs:.4f}, {ys:.4f}, {zs:.4f}) m\n")

print(f"Displacement Vector R = r_f - r_s")
print(f"  Rx = xf - xs = {xf} - {xs:.4f} = {Rx:.4f} m")
print(f"  Ry = yf - ys = {yf} - {ys:.4f} = {Ry:.4f} m")
print(f"  Rz = zf - zs = {zf} - {zs:.4f} = {Rz:.4f} m\n")

print(f"Magnitude of Displacement Vector R = sqrt(Rx² + Ry² + Rz²)")
print(f"  R = sqrt({Rx:.4f}² + {Ry:.4f}² + {Rz:.4f}²) = {R_mag:.4f} m\n")

print(f"Current Density Vector K in Cartesian coordinates (Kx, Ky, Kz)")
print(f"  Kx = -K_theta * sin(θ) = -{K_theta} * sin({theta_s_deg}°) = {Kx:.4f} A/m")
print(f"  Ky =  K_theta * cos(θ) =  {K_theta} * cos({theta_s_deg}°) = {Ky:.4f} A/m")
print(f"  Kz = 0 A/m\n")

print("Final Equation for Bz:")
print("Bz = (μ₀/4π) * ΔS * (Kx * Ry - Ky * Rx) / R³\n")
print("Substituting the values:")
print(f"Bz = ({mu_0_div_4pi:.1e}) * ({delta_S}) * (({Kx:.4f})*({Ry:.4f}) - ({Ky:.4f})*({Rx:.4f})) / ({R_mag:.4f})³")
print(f"Bz = ({mu_0_div_4pi:.1e}) * ({delta_S}) * ({Cz:.4f}) / ({R_mag**3:.4f})")
print(f"Bz = {Bz:.4e} Tesla\n")

print("--- Final Answer ---")
print(f"The magnitude of the longitudinal component of the magnetic field Bz is: {Bz:.3e} T")
<<<1.202e-05>>>