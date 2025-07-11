import math

# Given constants
P_GW = 1  # Luminosity in GW
S = 10  # Area of the cell in m^2
T = 12 * 3600  # Orbital period in seconds
M = 7.35e22  # Lunar mass in kg
R = 1738e3  # Lunar radius in m
G = 6.67e-11  # Gravitational constant in m^3 kg^-1 s^-2
A_mirror = 1e6 # Mirror area in m^2

# Convert P to Watts
P = P_GW * 1e9

# Step 1: Calculate the semi-major axis 'a' of the orbit using Kepler's Third Law
# T^2 = (4 * pi^2 / (G * M)) * a^3
# a^3 = (T^2 * G * M) / (4 * pi^2)
a_cubed = (T**2 * G * M) / (4 * math.pi**2)
a = a_cubed**(1/3)

# Step 2: Assume the distance between satellites d_XY is equal to the semi-major axis 'a'
d_XY = a

# Step 3: Calculate the power P' incident on the cell
# The formula derived from the virtual source model is P' = P * S / (4 * pi * d_XY^2)
P_prime = (P * S) / (4 * math.pi * d_XY**2)

# Step 4: Convert the result to microwatts and format the output
P_prime_microwatts = P_prime * 1e6

# Print the final equation with all the numerical values
print("Step 1: Calculate the semi-major axis 'a'.")
print(f"a = ( (T^2 * G * M) / (4 * pi^2) )^(1/3)")
print(f"a = ( ({T}^2 * {G} * {M}) / (4 * {math.pi:.5f}^2) )^(1/3)")
print(f"a = {a:.1f} m\n")

print("Step 2: Assume the distance between satellites equals the semi-major axis.")
print(f"d_XY = a = {d_XY:.1f} m\n")

print("Step 3: Calculate the final power P' using the simplified formula.")
print("P' = (P * S) / (4 * pi * d_XY^2)")
print(f"P' = ({P:.1e} W * {S} m^2) / (4 * {math.pi:.5f} * {d_XY:.1f}^2 m^2)")
print(f"P' = {P_prime:.4e} W\n")

print("Step 4: Convert to microwatts and round to one decimal place.")
print(f"P' (in microwatts) = {P_prime_microwatts:.1f} µW")

# Final answer in the required format
print(f"\nFinal Answer:")
print(f"The total power P' incident on the cell is {P_prime_microwatts:.1f} microwatts.")
print(f"<<<{P_prime_microwatts:.1f}>>>")