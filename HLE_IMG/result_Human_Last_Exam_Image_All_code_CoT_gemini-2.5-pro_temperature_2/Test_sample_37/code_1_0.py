import math

# --- Given Parameters ---
# Permeability of free space (mu_0)
mu_0 = 4 * math.pi * 1e-7  # T*m/A
# Azimuthal surface current density magnitude (K)
# Note: The problem states J = 40 A/m^2. Since the current is on a surface,
# this is interpreted as a surface current density K = 40 A/m.
# The calculation assumes a unit area source element (dA = 1 m^2) to find a finite B field.
K_mag = 40.0  # A/m
# Cylinder radius (r)
r_s = 0.34  # m
# Source point in cylindrical coordinates (r, theta, z)
theta_s_deg = 60.0  # degrees
z_s = 0.5  # m
# Target field point in Cartesian coordinates (x_f, y_f, z_f)
x_f = 0.1  # m
y_f = -0.1  # m
z_f = 0.2  # m

# --- Step 1: Coordinate Transformations ---
# Convert source angle to radians
theta_s_rad = math.radians(theta_s_deg)

# Convert source point to Cartesian coordinates
x_s = r_s * math.cos(theta_s_rad)
y_s = r_s * math.sin(theta_s_rad)
# z_s is already in Cartesian

# Convert azimuthal current vector K to Cartesian coordinates
# K_vec = K_mag * (-sin(theta) i + cos(theta) j)
Kx = -K_mag * math.sin(theta_s_rad)
Ky = K_mag * math.cos(theta_s_rad)

# --- Step 2: Calculate Displacement Vector R and its Magnitude ---
# Displacement vector R = r_f - r_s
Rx = x_f - x_s
Ry = y_f - y_s
Rz = z_f - z_s

# Magnitude of R
R_mag_sq = Rx**2 + Ry**2 + Rz**2
R_mag = math.sqrt(R_mag_sq)
R_mag_cubed = R_mag**3

# --- Step 3: Apply Biot-Savart Law for the z-component ---
# The z-component of K x R is (Kx*Ry - Ky*Rx)
cross_product_z = (Kx * Ry) - (Ky * Rx)

# Biot-Savart constant
biot_savart_const = mu_0 / (4 * math.pi)

# Calculate the longitudinal component of the magnetic field (Bz)
Bz = biot_savart_const * cross_product_z / R_mag_cubed

# --- Step 4: Output the results ---
print("--- Calculation Steps ---")
print(f"Source Point (Cartesian): (x_s, y_s, z_s) = ({x_s:.4f}, {y_s:.4f}, {z_s:.4f}) m")
print(f"Current Vector (Cartesian): (Kx, Ky) = ({Kx:.4f}, {Ky:.4f}) A/m")
print(f"Displacement Vector R: (Rx, Ry, Rz) = ({Rx:.4f}, {Ry:.4f}, {Rz:.4f}) m")
print(f"Magnitude of R: {R_mag:.4f} m")
print(f"Magnitude of R cubed: {R_mag_cubed:.4f} m^3")
print("\n--- Final Equation ---")
# Remember to print the full equation as requested
# The formula for Bz is (mu_0 / 4pi) * (Kx*Ry - Ky*Rx) / R^3
print(f"B_z = (1e-7) * (({Kx:.4f}) * ({Ry:.4f}) - ({Ky:.4f}) * ({Rx:.4f})) / {R_mag_cubed:.4f}")
print(f"B_z = (1e-7) * ({cross_product_z:.4f}) / {R_mag_cubed:.4f}")
print("\n--- Final Answer ---")
print(f"The magnitude of the longitudinal component of the magnetic field B_z is: {abs(Bz):.4e} T")

# For automated checking
final_answer = abs(Bz)