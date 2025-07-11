import numpy as np

# --- 1. Define constants and given values ---

# Permeability of free space constant (μ₀/4π) in H/m
mu0_4pi = 1e-7

# Target field point in Cartesian coordinates (xf, yf, zf) in meters
r_f = np.array([0.1, -0.1, 0.2])

# Source point in Cylindrical coordinates (r, θ, z)
r_s_cyl = np.array([0.34, 60.0, 0.5])  # θ in degrees

# Azimuthal surface current density K_theta in A/m (assuming A/m^2 was a typo)
K_theta = 40.0

# --- 2. Convert source point to Cartesian coordinates ---

# Convert angle from degrees to radians
theta_rad = np.deg2rad(r_s_cyl[1])

# Calculate Cartesian coordinates (xs, ys, zs)
x_s = r_s_cyl[0] * np.cos(theta_rad)
y_s = r_s_cyl[0] * np.sin(theta_rad)
z_s = r_s_cyl[2]
r_s = np.array([x_s, y_s, z_s])

# --- 3. Calculate displacement vector R and its magnitude ---

# Displacement vector R = r_f - r_s
R_vec = r_f - r_s

# Magnitude of R
R_mag = np.linalg.norm(R_vec)

# --- 4. Define the current density vector K in Cartesian coordinates ---

# The unit vector in the azimuthal direction (θ_hat) in Cartesian coordinates is (-sinθ, cosθ, 0)
theta_hat = np.array([-np.sin(theta_rad), np.cos(theta_rad), 0])

# Current density vector K
K_vec = K_theta * theta_hat

# --- 5. Calculate the cross product K x R ---

cross_product = np.cross(K_vec, R_vec)

# --- 6. Calculate the longitudinal magnetic field component Bz ---

# Assuming a unit surface area element dS = 1 m^2 to resolve the ill-posed problem
# The z-component of the magnetic field B is (μ₀/4π) * (K x R)_z / R^3
Bz = mu0_4pi * cross_product[2] / (R_mag**3)

# --- 7. Print the results step-by-step ---

print("Step-by-Step Calculation for Bz:\n")
print(f"1. Source point (Cartesian): r_s = ({r_s[0]:.4f}, {r_s[1]:.4f}, {r_s[2]:.4f}) m")
print(f"2. Field point (Cartesian):  r_f = ({r_f[0]:.4f}, {r_f[1]:.4f}, {r_f[2]:.4f}) m")
print(f"3. Displacement vector: R = r_f - r_s = ({R_vec[0]:.4f}, {R_vec[1]:.4f}, {R_vec[2]:.4f}) m")
print(f"4. Magnitude of R: |R| = {R_mag:.4f} m")
print(f"5. Current density vector (Cartesian): K = ({K_vec[0]:.4f}, {K_vec[1]:.4f}, {K_vec[2]:.4f}) A/m")

print("\nFormula for the z-component of the magnetic field:")
print("Bz = (μ₀/4π) * (Kx * Ry - Ky * Rx) / |R|³")

print("\nPlugging in the numbers:")
# The format instruction requires printing each number in the equation.
print(f"Bz = ({mu0_4pi:.1e}) * (({K_vec[0]:.4f}) * ({R_vec[1]:.4f}) - ({K_vec[1]:.4f}) * ({R_vec[0]:.4f})) / ({R_mag:.4f})³")
term1 = K_vec[0] * R_vec[1]
term2 = K_vec[1] * R_vec[0]
print(f"Bz = ({mu0_4pi:.1e}) * ({term1:.4f} - ({term2:.4f})) / ({R_mag**3:.6f})")
print(f"Bz = ({mu0_4pi:.1e}) * ({cross_product[2]:.4f}) / ({R_mag**3:.6f})")

print("\nFinal Result:")
print(f"The magnitude of the longitudinal component of the magnetic field Bz is: {abs(Bz):.4e} T")
print(f"\n<<< {abs(Bz):.4e} >>>")
