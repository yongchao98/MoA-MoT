import math

# Given constants in SI units
P = 1.0 * 10**9  # Luminosity of the source in Watts (1 GW)
S = 10.0  # Area of the photovoltaic cell in m^2
M = 7.35 * 10**22  # Mass of the Moon in kg
R = 1738.0 * 10**3  # Radius of the Moon in meters (1738 km)
G = 6.67 * 10**-11  # Gravitational constant in m^3 kg^-1 s^-2
T = 12.0 * 3600.0  # Orbital period in seconds (12 hours)
A_mirror = 1.0 * 10**6 # Area of satellite mirrors in m^2 (1 km^2)

# --- Step 1: Orbital Mechanics ---

# Calculate the semi-major axis 'a' of the orbit using Kepler's Third Law
# T^2 = (4 * pi^2 * a^3) / (G * M)
# a^3 = (G * M * T^2) / (4 * pi^2)
a_cubed = (G * M * T**2) / (4 * math.pi**2)
a = a_cubed**(1/3)

# We assume the relevant altitude is the one for a circular orbit with radius a
# This is a simplifying assumption to make the problem solvable without the orbit's eccentricity
h = a - R

# --- Step 2: Power Calculation ---

# The power calculation is based on a simplified model for chained reflections.
# The intensity (power per unit area) in the beam at the end of a leg of length r_leg
# is given by I = P_source / (2 * pi * r_leg^2).
# This relation holds for each leg of the A -> X -> Y -> B path.
# Therefore, the intensity at point B is determined by the length of the last leg, r_YB.
# We assume the satellite Y is at the calculated altitude h, so r_YB = h.
#
# Intensity at B: I_B = P / (2 * pi * h^2)
# Power on cell: P_prime = I_B * S

# The final power P' incident on the cell is:
P_prime = (P * S) / (2 * math.pi * h**2)

# --- Step 3: Final Result ---

# Convert the power from Watts to microwatts
P_prime_microwatts = P_prime * 10**6

# Print the detailed calculation steps and the final answer
print("--- Calculation Steps ---")
print(f"1. Gravitational Parameter (GM) = {G*M:.4e} m^3/s^2")
print(f"2. Orbital Period (T) = {T} s")
print(f"3. Calculated semi-major axis (a) = {a/1000:.2f} km")
print(f"4. Satellite altitude (h = a - R) = {h/1000:.2f} km")
print(f"5. Intensity at B (I_B = P / (2 * pi * h^2)) = ({P:.1e} W) / (2 * {math.pi:.2f} * ({h:.2e} m)^2) = {P / (2 * math.pi * h**2):.4e} W/m^2")
print(f"6. Final power on cell (P' = I_B * S) = ({P / (2 * math.pi * h**2):.4e} W/m^2) * {S} m^2 = {P_prime:.4e} W")
print("\n--- Final Answer ---")
print(f"The total power P' incident on the cell is {P_prime_microwatts:.1f} microwatts.")
