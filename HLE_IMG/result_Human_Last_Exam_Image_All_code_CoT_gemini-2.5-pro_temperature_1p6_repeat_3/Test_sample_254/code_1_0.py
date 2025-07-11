import math

# Constants
P_source = 1e9  # W (1 GW)
Area_mirror = 1e6  # m^2 (1 km^2)
S_cell = 10.0  # m^2
T_orbit = 12 * 3600  # s
M_moon = 7.35e22  # kg
R_moon = 1738e3  # m
G = 6.67e-11  # m^3 kg^-1 s^-2

# Step 1: Calculate the orbital radius (semi-major axis 'a')
# T^2 / a^3 = 4 * pi^2 / (G * M) => a = (T^2 * G * M / (4 * pi^2))^(1/3)
pi = math.pi
a_cubed = (T_orbit**2 * G * M_moon) / (4 * pi**2)
a = a_cubed**(1/3)

print("Step 1: Orbital Parameters")
print(f"The orbital period is T = {T_orbit} s")
print(f"The calculated semi-major axis (orbital radius) is a = {a:.1f} m")
print("-" * 30)

# Step 2 & 3: Geometry and Angles
# Based on the reasoning, we assume a circular orbit and a 90-degree separation
# between the satellites (angle XOY = 90 deg).
# This makes triangle OXY a right isosceles triangle.
# The angle of incidence for both reflections is 22.5 degrees.
theta_i_rad = math.radians(22.5)

# Calculate the effective area of the mirrors
A_eff = Area_mirror * math.cos(theta_i_rad)

# Calculate the relevant distances
dist_AX = a - R_moon
dist_YB = dist_AX
# dist_XY is the hypotenuse of the isosceles right triangle OXY
dist_XY = a * math.sqrt(2)

print("Step 2 & 3: System Geometry")
print("Assuming a circular orbit and 90-degree satellite separation:")
print(f"The distance from source A to satellite X is AX = {dist_AX:.1f} m")
print(f"The distance from satellite X to satellite Y is XY = {dist_XY:.1f} m")
print(f"The distance from satellite Y to detector B is YB = {dist_YB:.1f} m")
print(f"The angle of incidence at each mirror is 22.5 degrees, so cos(theta_i) = {math.cos(theta_i_rad):.4f}")
print("-" * 30)

# Step 4: Calculate power transfer
# a. Power intercepted by Mirror X
P_X = (P_source * A_eff) / (4 * pi * dist_AX**2)
print("Step 4: Power Transfer Calculation")
print("a. Power on Mirror X:")
print(f"P_X = ({P_source:.1e} W * {A_eff:.1f} m^2) / (4 * pi * ({dist_AX:.1f} m)^2) = {P_X:.3f} W")

# b. Power intercepted by Mirror Y
# The beam from A reflected at X diverges. Its solid angle is (A_eff / AX^2).
# Area of this beam at Y is Area_beam_Y = (A_eff / AX^2) * XY^2
# Power on Y is P_Y = P_X * (Area_eff / Area_beam_Y) = P_X * (AX / XY)^2
P_Y = P_X * (dist_AX / dist_XY)**2
print("\nb. Power on Mirror Y:")
print(f"P_Y = P_X * (AX / XY)^2 = {P_X:.3f} W * ({dist_AX:.1f} m / {dist_XY:.1f} m)^2 = {P_Y:.3f} W")


# c. Power intercepted by Detector B
# The beam from Y diverges as if from virtual source A'.
# We need the distance from A' to Y, |A'Y|.
# By the law of cosines in triangle A'XY, where angle A'XY = 135 deg:
# |A'Y|^2 = |AX|^2 + |XY|^2 - 2|AX||XY|cos(135)
dist_A_prime_Y_sq = dist_AX**2 + dist_XY**2 - 2 * dist_AX * dist_XY * math.cos(math.radians(135))
dist_A_prime_Y = math.sqrt(dist_A_prime_Y_sq)

# Solid angle of the beam reflected from Y is (A_eff / |A'Y|^2)
# Area of this beam at B is Area_beam_B = (A_eff / |A'Y|^2) * |YB|^2
Area_beam_B = (A_eff / dist_A_prime_Y_sq) * dist_YB**2

# Power on the cell is P' = P_Y * (S_cell / Area_beam_B), since S_cell is smaller than Area_beam_B
P_prime = P_Y * (S_cell / Area_beam_B)

print("\nc. Power on Detector B:")
print(f"The distance from virtual source A' to Y is |A'Y| = {dist_A_prime_Y:.1f} m")
print(f"The area of the light beam at B is Area_beam_B = {Area_beam_B:.1f} m^2")
print(f"P' = P_Y * (S_cell / Area_beam_B) = {P_Y:.3f} W * ({S_cell:.1f} m^2 / {Area_beam_B:.1f} m^2) = {P_prime:.8f} W")

# Step 5: Final Answer
# Convert power to microwatts
P_prime_microwatts = P_prime * 1e6
print("-" * 30)
print("Step 5: Final Result")
print(f"The final power P' on the cell is {P_prime:.8f} W.")
print(f"In microwatts, P' = {P_prime_microwatts:.1f} µW.")
