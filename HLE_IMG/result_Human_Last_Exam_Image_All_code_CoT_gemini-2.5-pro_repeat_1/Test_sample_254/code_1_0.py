import math

# --- Given Data ---
# Luminosity of the source
P = 1e9  # Watts (1 GW)
# Area of the photovoltaic cell
S = 10.0  # m^2
# Area of each satellite's mirror
A_mirror = 1.0 * (1000**2)  # km^2 to m^2
# Orbital period
T = 12.0 * 3600.0  # hours to seconds
# Lunar mass
M = 7.35e22  # kg
# Lunar radius
R = 1738.0 * 1000.0  # km to m
# Gravitational constant
G = 6.67e-11  # m^3 kg^-1 s^-2

# --- Step 1 & 2: Calculate the orbital radius (r) ---
# T^2 = (4 * pi^2 * r^3) / (G * M)
# r^3 = (G * M * T^2) / (4 * pi^2)
pi = math.pi
r_cubed = (G * M * T**2) / (4 * pi**2)
r = r_cubed**(1/3)

print("--- Calculation Steps ---")
print(f"Gravitational Constant (G): {G:.2e} m^3 kg^-1 s^-2")
print(f"Lunar Mass (M): {M:.2e} kg")
print(f"Orbital Period (T): {T:.0f} s")
print(f"Calculated orbital radius (r): {r/1000:.1f} km")

# --- Step 3: Determine the geometry and distances for the virtual source model ---
# Based on the analysis, we use a circular orbit with antipodal satellites.
# Let the Moon's center be the origin (0).
# A is at R, B is at -R.
# X is at r, Y is at -r.
# 1. Virtual source A' from reflection at X:
#    Mirror is at x_m1 = r. Source is at x_s1 = R.
#    Position of A': x_v1 = 2*x_m1 - x_s1 = 2*r - R.
# 2. Virtual source A'' from reflection at Y:
#    Mirror is at x_m2 = -r. Source is virtual source A' at x_s2 = x_v1.
#    Position of A'': x_v2 = 2*x_m2 - x_s2 = 2*(-r) - (2*r - R) = -4*r + R.
# 3. Distance from A'' to B:
#    Detector B is at x_b = -R.
#    d_A_double_prime_B = abs(x_b - x_v2) = abs(-R - (-4*r + R)) = abs(4*r - 2*R)
d_A_double_prime_B = 4*r - 2*R

print(f"Lunar Radius (R): {R/1000:.1f} km")
print(f"Position of virtual source A': { (2*r-R)/1000:.1f} km")
print(f"Position of virtual source A'': {(-4*r+R)/1000:.1f} km")
print(f"Distance from final virtual source A'' to detector B: {d_A_double_prime_B/1000:.1f} km")


# --- Step 4 & 5: Calculate the final power P' ---
# The radiant intensity J of the source is P / (2*pi) because it radiates into a hemisphere.
# This radiant intensity is conserved for the virtual sources (assuming perfect reflection).
J = P / (2 * pi)
# The power P' on the detector is the flux from A'' at B, multiplied by the cell's area S.
# Flux = J / d^2
P_prime = J * S / (d_A_double_prime_B**2)

# Convert to microwatts
P_prime_microwatts = P_prime * 1e6

print("\n--- Final Calculation ---")
print(f"Source Luminosity (P): {P:.1e} W")
print(f"Source Radiant Intensity (J = P / 2π): {J:.2e} W/sr")
print(f"Detector Area (S): {S:.1f} m^2")
print(f"Final calculated power incident on the cell (P'): {P_prime:.3e} W")
print(f"Final power in microwatts: {P_prime_microwatts:.1f} μW")

# Final answer formatted as requested
# The final equation is P' = (P / (2 * π)) * S / (4r - 2R)²
# We print the components of this equation again for clarity.
print("\n--- Final Equation P' = (P / (2 * π)) * S / (4r - 2R)² ---")
print(f"P = {P:.1e} W")
print(f"S = {S:.1f} m²")
print(f"r = {r:.3e} m")
print(f"R = {R:.3e} m")
print(f"P' = ({P:.1e} / (2 * {pi:.4f})) * {S:.1f} / (4*{r:.3e} - 2*{R:.3e})²")
print(f"P' = {P_prime_microwatts:.1f} μW")

# Output the final numerical answer
# <<<3.6>>>