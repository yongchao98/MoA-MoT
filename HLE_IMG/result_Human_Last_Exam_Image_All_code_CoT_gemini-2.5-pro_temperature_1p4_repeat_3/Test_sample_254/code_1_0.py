import math

# Given constants
P_source = 1e9  # Luminosity of the source in Watts (1 GW)
S_cell = 10.0  # Area of the photovoltaic cell in m^2
M_moon = 7.35e22  # Mass of the Moon in kg
R_moon = 1738e3  # Radius of the Moon in meters (1738 km)
G = 6.67e-11  # Gravitational constant in m^3 kg^-1 s^-2
T_orbit = 12 * 3600  # Orbital period in seconds (12 hours)

# Step 1: Calculate the semi-major axis 'a' of the orbit using Kepler's Third Law
# T^2 / a^3 = 4 * pi^2 / (G * M) => a = (G * M * T^2 / (4 * pi^2))^(1/3)
# We assume a circular orbit for a unique solution, so the orbital radius r is equal to a.
a_cubed = (G * M_moon * T_orbit**2) / (4 * math.pi**2)
a = a_cubed**(1/3)
r_orbit = a

print(f"Calculations based on a simplified symmetric model (circular orbit):")
print(f"1. Orbital Parameters:")
print(f"  - Given Orbital Period T = {T_orbit} s")
print(f"  - Calculated Semi-major Axis a = {r_orbit / 1000:.1f} km")

# Step 2: Calculate the distance from the source A to satellite X
# Satellite X is at the zenith of A, so the distance is its altitude.
d_AX = r_orbit - R_moon
print(f"2. Light Path Distances:")
print(f"  - Altitude of satellite (distance A to X), d_AX = {d_AX / 1000:.1f} km")

# Step 3: Calculate the intensity of light at satellite X
# The isotropic source P radiates over a solid angle of 4*pi.
# The surface of the moon blocks half of this.
# Intensity I = Power / Area = P / (4 * pi * d^2)
intensity_at_X = P_source / (4 * math.pi * d_AX**2)
print(f"3. Power Calculation:")
print(f"  - Intensity from source at Satellite X, I_X = {intensity_at_X:.4g} W/m^2")

# Step 4: Calculate the power incident on the cell at B
# Neglecting diffraction and assuming perfect mirrors, the beam from Y to B has the same
# intensity as the light arriving at X.
# The cell at B lies horizontally, and in our symmetric model, the beam is incident normally.
# Power_on_cell P' = Intensity * Cell_Area
power_on_cell_W = intensity_at_X * S_cell
print(f"  - The intensity is conserved through the reflection chain (I_B = I_X).")
print(f"  - Power on cell P' = I_B * S = {intensity_at_X:.4g} W/m^2 * {S_cell} m^2 = {power_on_cell_W:.4g} W")

# Step 5: Convert the result to microwatts
power_on_cell_uW = power_on_cell_W * 1e6
print("\nFinal Result:")
print(f"The final power P' incident on the cell is {power_on_cell_uW:.1f} microwatts.")

# The final equation is P' = (P * S) / (4 * pi * (a - R)^2)
print("\nFinal Equation used:")
print(f"P' = ({P_source:.1e} W * {S_cell} m^2) / (4 * pi * ({r_orbit:.4g} m - {R_moon:.4g} m)^2)")
print(f"P' = ({P_source*S_cell:.1e}) / (4 * pi * ({d_AX:.4g})^2) = {power_on_cell_uW:.1f} uW")

<<<41.0>>>