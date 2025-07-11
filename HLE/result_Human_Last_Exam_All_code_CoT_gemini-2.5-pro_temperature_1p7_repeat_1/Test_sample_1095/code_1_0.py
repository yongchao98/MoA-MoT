import numpy as np

print("This script demonstrates why the condition k_r ∝ √l is necessary for rigid rotation in Bessel-Gauss (BG) beams.")
print("We calculate the rotation rate of a wave packet formed by superposing two BG modes (l1 and l2).")
print("Rigid rotation requires this rate to be constant, which is achieved if k_r(l)² is linear in l.")

# --- Define physical parameters for the example ---
wavelength = 633e-9  # Wavelength of a HeNe laser in meters
l1 = 3  # Topological charge of the first mode
l2 = 5  # Topological charge of the second mode
# Proportionality constant C for the condition k_r = C * √l
C = 1000.0  # Units of m^(-1/2)

# --- Calculations ---
# Total wavevector k
k = 2 * np.pi / wavelength

# Calculate the radial wavevectors k_r based on the k_r ∝ √l condition
k_r_l1 = C * np.sqrt(l1)
k_r_l2 = C * np.sqrt(l2)

# Calculate propagation constants k_z using the paraxial approximation: k_z ≈ k - k_r² / (2*k)
k_z_l1 = k - (k_r_l1**2) / (2 * k)
k_z_l2 = k - (k_r_l2**2) / (2 * k)

# The rotation rate of the interference pattern is dφ/dz = (k_z(l1) - k_z(l2)) / (l1 - l2)
rotation_rate = (k_z_l1 - k_z_l2) / (l1 - l2)

# Theoretical rate: From the formulas, the rate simplifies to -C² / (2*k)
theoretical_rate = -(C**2) / (2 * k)

# --- Print the results and equations ---
print("\n--- Example Calculation ---")
print(f"Parameters: λ = {wavelength:.2e} m, C = {C}, l1 = {l1}, l2 = {l2}")
print(f"Total wavevector k = 2π / {wavelength:.2e} m = {k:.4e} rad/m")

print("\n1. Radial wavevectors from k_r = C * √l:")
print(f"   k_r(l={l1}) = {C:.1f} * √{l1} = {k_r_l1:.4e} rad/m")
print(f"   k_r(l={l2}) = {C:.1f} * √{l2} = {k_r_l2:.4e} rad/m")

print("\n2. Propagation constants from k_z = k - k_r² / (2*k):")
print(f"   k_z(l={l1}) = {k:.4e} - ({k_r_l1:.4e})² / (2 * {k:.4e}) = {k_z_l1:.4e} rad/m")
print(f"   k_z(l={l2}) = {k:.4e} - ({k_r_l2:.4e})² / (2 * {k:.4e}) = {k_z_l2:.4e} rad/m")

print("\n3. Final rotation rate of the wave packet:")
print(f"   Rate = (k_z(l1) - k_z(l2)) / (l1 - l2)")
print(f"   Rate = ({k_z_l1:.4e} rad/m - {k_z_l2:.4e} rad/m) / ({l1} - {l2})")
print(f"   Rate = {rotation_rate:.4f} rad/m")

print(f"\nThis calculated rate matches the theoretical constant rate of -C²/(2*k) = {theoretical_rate:.4f} rad/m.")
print("Since the rate is constant, the condition k_r ∝ √l ensures rigid rotation.")
