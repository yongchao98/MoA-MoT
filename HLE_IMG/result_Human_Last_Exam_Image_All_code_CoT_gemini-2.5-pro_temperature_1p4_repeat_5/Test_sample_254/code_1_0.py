import math

# Given constants
P_source = 1e9  # Power of the source in Watts (1 GW)
S_cell = 10.0   # Area of the photovoltaic cell in m^2
T_period = 12 * 3600  # Orbital period in seconds (12 hours)
G = 6.67e-11  # Gravitational constant in kg^-1 m^3 s^-2
M_moon = 7.35e22  # Mass of the Moon in kg
R_moon = 1738 * 1000  # Radius of the Moon in meters (1738 km)
Area_mirror = 1e6 # Area of the mirror in m^2 (1 km^2)

# Step 1: Calculate the semi-major axis 'a' using Kepler's Third Law
# T^2 = (4 * pi^2 / (G * M)) * a^3
a_cubed = (T_period**2 * G * M_moon) / (4 * math.pi**2)
a = a_cubed**(1/3.0)

# Step 2: Determine the total unfolded path length 'L'
# The geometry implies a circular orbit (r=a) and d_XY = a
d_AX = a - R_moon
d_XY = a
d_YB = a - R_moon
L_total = d_AX + d_XY + d_YB
# Which simplifies to L = 3*a - 2*R

# Step 3: Calculate the final power P' incident on the cell
# P' = P_source * S_cell / (4 * pi * L_total^2)
P_prime = (P_source * S_cell) / (4 * math.pi * L_total**2)

# Convert final power to microwatts
P_prime_microwatts = P_prime * 1e6

# Print the results step-by-step
print("--- Calculation Steps ---")
print(f"1. Orbital period T = {T_period:.0f} s")
print(f"2. Gravitational constant G = {G:.2e} kg^-1 m^3 s^-2")
print(f"3. Mass of the Moon M = {M_moon:.2e} kg")
print(f"4. Radius of the Moon R = {R_moon:.0f} m")

# Print the formula and calculation for the semi-major axis 'a'
print("\nCalculating the orbital radius 'a' (semi-major axis for a circular orbit):")
print(f"a = (T^2 * G * M / (4 * pi^2))^(1/3)")
print(f"a = ({T_period:.0f}^2 * {G:.2e} * {M_moon:.2e} / (4 * pi^2))^(1/3)")
print(f"a = {a:.3e} m")

# Print the formula and calculation for the total path length 'L'
print("\nCalculating the total unfolded light path length 'L':")
print(f"L = 3*a - 2*R")
print(f"L = 3 * {a:.3e} - 2 * {R_moon:.0f}")
print(f"L = {L_total:.3e} m")

# Print the final formula and calculation for the power P'
print("\nCalculating the final power P' on the cell:")
print(f"P' = P_source * S_cell / (4 * pi * L^2)")
print(f"P' = {P_source:.0e} W * {S_cell:.1f} m^2 / (4 * pi * {L_total:.3e}^2 m^2)")
print(f"P' = {P_prime:.3e} W")

# Print the final answer in the required format
print("\n--- Final Answer ---")
print(f"The total power incident on the cell is {P_prime_microwatts:.1f} microwatts.")
<<<3.6>>>