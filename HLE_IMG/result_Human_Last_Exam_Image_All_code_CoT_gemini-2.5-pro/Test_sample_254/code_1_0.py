import math

# --- Given Data (in SI units) ---
P = 1e9  # Power of the source in Watts (1 GW)
S = 10.0  # Area of the photovoltaic cell in m^2
M = 7.35e22  # Mass of the Moon in kg
R = 1738e3  # Radius of the Moon in meters (1738 km)
G = 6.67e-11  # Gravitational constant in kg^-1 m^3 s^-2
T = 12 * 3600  # Orbital period in seconds (12 hours)

# --- Step 1: Calculate the semi-major axis 'a' of the orbit ---
# Using Kepler's Third Law: T^2 = (4 * pi^2 / (G * M)) * a^3
# Rearranging for a: a = (G * M * T^2 / (4 * pi^2))^(1/3)
four_pi_squared = 4 * math.pi**2
a_cubed = (G * M * T**2) / four_pi_squared
a = a_cubed**(1/3)

# --- Step 2: Calculate the total unfolded path length d_total ---
# d_total = 4a - 2R
d_total = 4 * a - 2 * R

# --- Step 3: Calculate the final power P' on the cell ---
# The intensity at distance d_total is I = P / (2 * pi * d_total^2)
# The power on the cell is P' = I * S
# P' = (P * S) / (2 * pi * d_total^2)
P_prime_watts = (P * S) / (2 * math.pi * d_total**2)

# --- Step 4: Convert the result to microwatts ---
P_prime_microwatts = P_prime_watts * 1e6

# --- Output the results ---
print("This script calculates the power received by a photovoltaic cell on the Moon.")
print("The calculation follows a specific geometric interpretation of the problem.")
print("\n--- Final Equation ---")
print("P' = (P * S) / (2 * pi * (4*a - 2*R)^2)")

print("\n--- Values Used ---")
print(f"Source Power (P): {P:.1e} W")
print(f"Cell Area (S): {S:.1f} m^2")
print(f"Calculated Semi-Major Axis (a): {a:.3e} m")
print(f"Lunar Radius (R): {R:.3e} m")

print("\n--- Calculation Steps ---")
print(f"Total unfolded path length (d_total = 4a - 2R): {d_total:.3e} m")
print(f"Power on cell (P'): {P_prime_watts:.3e} W")

print("\n--- Final Answer ---")
print(f"The total power P' incident on the cell is: {P_prime_microwatts:.1f} microwatts")

# The final numerical answer as requested
# print(f"<<<{P_prime_microwatts:.1f}>>>")