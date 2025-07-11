import math

# Step 1: Define Constants and Given Values
mu_0 = 4 * math.pi * 1e-7  # Permeability of free space in H/m or T*m/A
mu0_over_4pi = 1e-7       # The constant factor in Biot-Savart Law

# Source point in cylindrical coordinates (r, theta, z)
r_s_cyl = 0.34  # m
theta_s_deg = 60.0  # degrees
z_s = 0.5  # m

# Target field point in Cartesian coordinates (x, y, z)
xf, yf, zf = 0.1, -0.1, 0.2

# Azimuthal surface current density K (assuming A/m from a typo in A/m^2)
K_theta_mag = 40.0  # A/m

# --- Calculations ---

# Step 2: Convert source point and current vector to Cartesian coordinates
# Convert theta from degrees to radians
theta_s_rad = math.radians(theta_s_deg)

# Source point in Cartesian coordinates
xs = r_s_cyl * math.cos(theta_s_rad)
ys = r_s_cyl * math.sin(theta_s_rad)
# zs is already given

# The azimuthal unit vector theta_hat in Cartesian coordinates is (-sin(theta), cos(theta), 0)
Kx = -K_theta_mag * math.sin(theta_s_rad)
Ky = K_theta_mag * math.cos(theta_s_rad)
Kz = 0.0

# Step 3: Calculate the displacement vector R = r_f - r_s
Rx = xf - xs
Ry = yf - ys
Rz = zf - z_s

# Calculate the magnitude of R
R_mag_sq = Rx**2 + Ry**2 + Rz**2
R_mag = math.sqrt(R_mag_sq)
R_mag_cubed = R_mag**3

# Step 4 & 5: Calculate the z-component of the cross product (K x R)
# (K x R)_z = Kx * Ry - Ky * Rx
cross_product_z = Kx * Ry - Ky * Rx

# Step 6: Calculate the longitudinal component of the magnetic field (Bz)
Bz = mu0_over_4pi * cross_product_z / R_mag_cubed

# Step 7: Output the results and the final equation
print("Calculation of the longitudinal magnetic field component (Bz):")
print("----------------------------------------------------------------")
print(f"Source Point (Cartesian): xs={xs:.4f} m, ys={ys:.4f} m, zs={zs:.4f} m")
print(f"Target Point (Cartesian): xf={xf:.4f} m, yf={yf:.4f} m, zf={zf:.4f} m")
print(f"Current Density Vector (Cartesian): Kx={Kx:.4f} A/m, Ky={Ky:.4f} A/m, Kz={Kz:.4f} A/m")
print(f"Displacement Vector R: Rx={Rx:.4f} m, Ry={Ry:.4f} m, Rz={Rz:.4f} m")
print("----------------------------------------------------------------")
print("The longitudinal component of the magnetic field is given by:")
print("Bz = (μ₀/4π) * (K × R)_z / |R|³")
print("\nPlugging in the calculated values:")
print(f"Bz = ({mu0_over_4pi:.1e} T·m/A) * ({cross_product_z:.4f} A·m) / ({R_mag_cubed:.4f} m³)")

print(f"\nFinal Result:")
print(f"Bz = {Bz:.4e} T")
print("----------------------------------------------------------------")

<<<1.2017e-05>>>