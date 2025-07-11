import math

# --- Step 1: Define constants and parameters ---
# Permeability of free space constant
mu0_4pi = 1e-7  # mu_0 / (4 * pi) in T*m/A

# Target field point in Cartesian coordinates (m)
xf, yf, zf = 0.1, -0.1, 0.2

# Source point in Cylindrical coordinates (m, degrees)
r_s = 0.34
theta_s_deg = 60.0
z_s = 0.5

# Azimuthal current term. Interpreted as the magnitude of a current moment vector M
# to ensure dimensional consistency. Assumed units: A*m^2.
M_mag = 40.0

# --- Step 2: Convert source point to Cartesian coordinates ---
# Convert theta from degrees to radians
theta_s_rad = math.radians(theta_s_deg)

# Calculate Cartesian coordinates of the source point
xs = r_s * math.cos(theta_s_rad)
ys = r_s * math.sin(theta_s_rad)

# --- Step 3: Calculate the displacement vector R ---
Rx = xf - xs
Ry = yf - ys
Rz = zf - zs

# --- Step 4: Calculate the magnitude of R ---
R_mag_sq = Rx**2 + Ry**2 + Rz**2
R_mag = math.sqrt(R_mag_sq)
R_mag_cubed = R_mag**3

# --- Step 5: Define the source current vector M in Cartesian coordinates ---
# The direction of M is azimuthal (theta_hat)
# The Cartesian components of the azimuthal unit vector are (-sin(theta), cos(theta), 0)
Mx = -M_mag * math.sin(theta_s_rad)
My = M_mag * math.cos(theta_s_rad)
Mz = 0

# --- Step 6: Calculate the z-component of the cross product (M x R) ---
Cross_product_z = Mx * Ry - My * Rx

# --- Step 7: Apply the Biot-Savart Law to find Bz ---
Bz = mu0_4pi * Cross_product_z / R_mag_cubed

# --- Step 8: Print the calculation and results ---
print("--- Calculation of the Magnetic Field Component B_z ---")
print(f"Target point (xf, yf, zf) = ({xf}, {yf}, {zf}) m")
print(f"Source point (r, θ, z) = ({r_s}, {theta_s_deg}, {z_s}) m")
print(f"Source point in Cartesian (xs, ys, zs) = ({xs:.5f}, {ys:.5f}, {zs:.1f}) m")
print(f"Displacement vector R = (Rx, Ry, Rz) = ({Rx:.5f}, {Ry:.5f}, {Rz:.1f}) m")
print(f"Current vector M = (Mx, My, Mz) = ({Mx:.5f}, {My:.5f}, {Mz:.1f}) A·m²")

print("\n--- Final Equation ---")
print("B_z = (μ₀/4π) * (M x R)_z / |R|³")
print(f"B_z = ({mu0_4pi}) T·m/A * ({Cross_product_z:.5f}) A·m³ / ({R_mag:.5f})³ m³")
print(f"B_z = ({mu0_4pi}) T·m/A * ({Cross_product_z:.5f}) A·m³ / ({R_mag_cubed:.5f}) m³")
print(f"B_z = {Bz:.5e} T")
<<<1.20172e-05>>>