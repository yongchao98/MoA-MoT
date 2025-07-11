import math

# --- Given Data ---
# Power of the source
P = 1e9  # GW to W
# Area of the photovoltaic cell
S = 10.0  # m^2
# Area of each satellite's mirror
A_mirror = 1e6 # 1 km^2 to m^2
# Orbital period
T = 12 * 3600  # hours to seconds
# Lunar mass
M = 7.35e22  # kg
# Lunar radius
R = 1738e3  # km to m
# Gravitational constant
G = 6.67e-11  # kg^-1 m^3 s^-2

# --- Step 1: Calculate the semi-major axis 'a' of the orbit ---
# Using Kepler's Third Law: T^2 = (4 * pi^2 / (G * M)) * a^3
# a = (G * M * T^2 / (4 * pi^2))^(1/3)
four_pi_sq = 4 * math.pi**2
a_cubed = (G * M * T**2) / four_pi_sq
a = a_cubed**(1/3)

# --- Step 2: Calculate the orbital altitude 'h' ---
# Assuming a circular orbit for simplification, where radius r = a
# h = a - R
h = a - R

# --- Step 3: Calculate the incident power P' on the cell ---
# The derived formula for power transfer is P' = (P * S) / (2 * pi * h^2)
P_prime = (P * S) / (2 * math.pi * h**2)

# --- Step 4: Convert to microwatts and format the output ---
P_prime_microwatts = P_prime * 1e6

# --- Print the results and the final equation with values ---
print("--- Calculation Steps ---")
print(f"1. Semi-major axis (a):")
print(f"   a = (G * M * T^2 / (4 * pi^2))^(1/3)")
print(f"   a = ({G:.2e} * {M:.2e} * {T:.2e}^2 / (4 * {math.pi:.4f}^2))^(1/3)")
print(f"   a = {a:.3e} m")
print("\n2. Orbital altitude (h):")
print(f"   h = a - R")
print(f"   h = {a:.3e} m - {R:.3e} m")
print(f"   h = {h:.3e} m")
print("\n3. Power on cell (P'):")
print(f"   P' = (P * S) / (2 * pi * h^2)")
print(f"   P' = ({P:.1e} W * {S:.1f} m^2) / (2 * {math.pi:.4f} * ({h:.3e} m)^2)")
print(f"   P' = {P_prime:.3e} W")
print("\n--- Final Answer ---")
print(f"The total power incident on the cell is {P_prime_microwatts:.1f} microwatts.")

# The final numerical answer in the required format
final_answer = f"{P_prime_microwatts:.1f}"
# print(f"<<<{final_answer}>>>")