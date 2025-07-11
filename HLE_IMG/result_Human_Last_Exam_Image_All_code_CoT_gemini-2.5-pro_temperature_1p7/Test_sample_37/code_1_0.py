import math

# --- Given Parameters ---
# Permeability of free space constant
mu_0_div_4pi = 1e-7  # H/m or T*m/A

# Source point in cylindrical coordinates
r_s_cyl = 0.34  # m
theta_s_deg = 60.0  # degrees
z_s = 0.5  # m

# Target field point in Cartesian coordinates
x_f = 0.1   # m
y_f = -0.1  # m
z_f = 0.2   # m

# Magnitude of the effective current element (interpreted from the problem statement)
Idl_mag = 40.0  # A*m

# --- Step 1: Convert source point to Cartesian coordinates ---
theta_s_rad = math.radians(theta_s_deg)
x_s = r_s_cyl * math.cos(theta_s_rad)
y_s = r_s_cyl * math.sin(theta_s_rad)
# z_s remains the same

# --- Step 2: Calculate the displacement vector R = r_f - r_s ---
R_x = x_f - x_s
R_y = y_f - y_s
R_z = z_f - z_s

# --- Step 3: Calculate the magnitude of R cubed (R^3) ---
R_squared = R_x**2 + R_y**2 + R_z**2
R = math.sqrt(R_squared)
R_cubed = R**3

# --- Step 4: Express the current element vector Idl in Cartesian coordinates ---
# The azimuthal unit vector in Cartesian is (-sin(theta), cos(theta), 0)
Idl_x = -Idl_mag * math.sin(theta_s_rad)
Idl_y = Idl_mag * math.cos(theta_s_rad)
# Idl_z = 0

# --- Step 5: Calculate the z-component of the cross product (Idl x R) ---
cross_product_z = Idl_x * R_y - Idl_y * R_x

# --- Step 6: Calculate the longitudinal component of the magnetic field B_z ---
B_z = mu_0_div_4pi * cross_product_z / R_cubed

# --- Final Output ---
print("The calculation for the longitudinal component of the magnetic field (B_z) is based on the Biot-Savart law:")
print("B_z = (μ₀/4π) * ((Idl × R)_z) / R³\n")
print(f"The equation with the calculated values is:")
print(f"B_z = ({mu_0_div_4pi}) * ({cross_product_z:.5f}) / ({R_cubed:.5f})")
print(f"B_z = {B_z:.4e} T")

# For clarity, here are the intermediate vector components:
print("\n--- Intermediate values ---")
print(f"Source point r_s = ({x_s:.3f}, {y_s:.3f}, {z_s:.3f}) m")
print(f"Target point r_f = ({x_f:.3f}, {y_f:.3f}, {z_f:.3f}) m")
print(f"Displacement vector R = ({R_x:.3f}, {R_y:.3f}, {R_z:.3f}) m")
print(f"Current element Idl = ({Idl_x:.3f}, {Idl_y:.3f}, 0.000) A·m")
<<<1.2017e-05>>>