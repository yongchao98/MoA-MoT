import math

# Define the given parameters
# mass of the body in kg
m = 0.20
# mass of the guide in kg
M = 0.80
# radius of the circular arcs in meters (20 cm = 0.2 m)
R = 0.20
# length of the straight section in meters (50 cm = 0.5 m)
d = 0.50

# The principle of conservation of the center of mass position is used.
# The horizontal displacement of the guide (delta_X_M) is given by the formula:
# delta_X_M = - (m / (m + M)) * (d + 2 * R)
# The negative sign indicates that the guide moves in the opposite direction
# to the horizontal movement of the mass m.

# Calculate the numerator of the mass ratio
mass_ratio_num = m

# Calculate the denominator of the mass ratio
mass_ratio_den = m + M

# Calculate the total horizontal distance the mass travels relative to the guide
relative_distance = d + 2 * R

# Calculate the horizontal displacement of the guide
delta_X_M = - (mass_ratio_num / mass_ratio_den) * relative_distance

print("The horizontal displacement of the guide (ΔX_M) is calculated using the conservation of the center of mass.")
print("The formula is: ΔX_M = - (m / (m + M)) * (d + 2*R)")
print("\nPlugging in the values:")
print(f"ΔX_M = - ({mass_ratio_num} / ({m} + {M})) * ({d} + 2*{R})")
print(f"ΔX_M = - ({mass_ratio_num} / {mass_ratio_den}) * ({relative_distance})")
print(f"ΔX_M = {delta_X_M:.2f} meters")
print("\nThis means the guide moves 18 cm in the direction opposite to the block's main horizontal movement.")

# The final numerical answer in meters
final_answer = delta_X_M
print(f"<<<{final_answer}>>>")