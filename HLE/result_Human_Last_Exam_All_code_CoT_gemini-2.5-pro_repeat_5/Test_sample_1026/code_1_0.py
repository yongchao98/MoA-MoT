# Define the problem parameters from the user's request
m = 0.20  # Mass of the body in kg
M = 0.80  # Mass of the guide in kg
d_cm = 50  # Length of the straight section in cm

# Convert the length to SI units (meters)
d = d_cm / 100.0

# The principle of conservation of the center of mass states that
# the displacement of the guide (Δx_M) can be found using the formula:
# Δx_M = - (m / (m + M)) * Δx_m_rel
# where Δx_m_rel is the horizontal displacement of the mass relative to the guide.

# In this problem, the relative horizontal displacement is the length of the straight section, d.
delta_x_m_relative = d

# Calculate the total mass of the system
total_mass = m + M

# Calculate the horizontal displacement of the guide
delta_x_M = - (m / total_mass) * delta_x_m_relative

print("This script calculates the horizontal displacement of the guide based on the conservation of the center of mass.")
print("\nThe final equation for the guide's displacement (Δx_guide) is:")
# The user requested to see each number in the final equation
print(f"Δx_guide = - (mass_body / (mass_body + mass_guide)) * relative_displacement_of_mass")
print(f"Δx_guide = - ({m} / ({m} + {M})) * {delta_x_m_relative}")

# Print the final calculated result
print(f"\nResult:")
print(f"The horizontal displacement of the guide is {delta_x_M:.2f} m.")
print("The negative sign indicates the guide moves to the left, opposite to the mass's general direction of travel.")
