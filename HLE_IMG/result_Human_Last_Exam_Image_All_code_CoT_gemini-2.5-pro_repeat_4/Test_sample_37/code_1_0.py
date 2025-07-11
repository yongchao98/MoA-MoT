import math

# --- Given Constants and Values ---
# Permeability of free space constant (mu_0 / 4pi) in T*m/A
mu0_over_4pi = 1e-7

# Target field point in Cartesian coordinates (m)
xf, yf, zf = 0.1, -0.1, 0.2

# Source point in Cylindrical coordinates (m, degrees)
r_s = 0.34
theta_s_deg = 60
z_s = 0.5

# Magnitude of the current element (A*m)
# Interpreted from the given "azimuthal current density" of 40 A/m^2
Idl_mag = 40.0

# --- Step 1: Convert source point to Cartesian coordinates ---
# Convert theta from degrees to radians
theta_s_rad = math.radians(theta_s_deg)

# Calculate Cartesian coordinates of the source point
xs = r_s * math.cos(theta_s_rad)
ys = r_s * math.sin(theta_s_rad)
# zs is the same as in cylindrical coordinates

# --- Step 2: Calculate the displacement vector R ---
Rx = xf - xs
Ry = yf - ys
Rz = zf - z_s
R_mag = math.sqrt(Rx**2 + Ry**2 + Rz**2)

# --- Step 3: Determine the current element vector Idl in Cartesian coordinates ---
# The direction is azimuthal (theta_hat).
# theta_hat = -sin(theta) * x_hat + cos(theta) * y_hat
Idl_x = Idl_mag * -math.sin(theta_s_rad)
Idl_y = Idl_mag * math.cos(theta_s_rad)
Idl_z = 0  # No z-component for azimuthal current

# --- Step 4 & 5: Calculate the z-component of the cross product and Bz ---
# (Idl x R)_z = Idl_x * Ry - Idl_y * Rx
cross_product_z = Idl_x * Ry - Idl_y * Rx

# Calculate the longitudinal (z-component) of the magnetic field
Bz = mu0_over_4pi * cross_product_z / (R_mag**3)

# --- Step 6: Print the results ---
print("Calculation of the longitudinal magnetic field component (Bz)")
print("-" * 60)
print(f"Target point (xf, yf, zf) = ({xf}, {yf}, {zf}) m")
print(f"Source point (r, theta, z) = ({r_s}, {theta_s_deg} deg, {z_s}) m")
print("-" * 60)
print("Intermediate Calculations:")
print(f"Source point in Cartesian (xs, ys, zs) = ({xs:.4f}, {ys:.4f}, {zs:.4f}) m")
print(f"Displacement vector R = (Rx, Ry, Rz) = ({Rx:.4f}, {Ry:.4f}, {Rz:.4f}) m")
print(f"Magnitude of R = {R_mag:.4f} m")
print(f"Current element Idl = (Idl_x, Idl_y, Idl_z) = ({Idl_x:.4f}, {Idl_y:.4f}, {Idl_z:.4f}) A.m")
print(f"z-component of cross product (Idl x R)z = {cross_product_z:.4f} A.m^2")
print("-" * 60)
print("Final Equation (Biot-Savart Law for Bz):")
print("Bz = (mu0 / 4pi) * (Idl_x * Ry - Idl_y * Rx) / |R|^3")
print("\nSubstituting the values:")
print(f"Bz = ({mu0_over_4pi}) * (({Idl_x:.4f}) * ({Ry:.4f}) - ({Idl_y:.4f}) * ({Rx:.4f})) / ({R_mag:.4f})^3")
print(f"Bz = ({mu0_over_4pi}) * ({cross_product_z:.4f}) / ({R_mag**3:.4f})")
print(f"Bz = {Bz:.4e} T")
print("-" * 60)
print(f"The magnitude of the longitudinal component of the magnetic field is {abs(Bz):.4e} Tesla.")

# The final answer in the required format
final_answer = abs(Bz)