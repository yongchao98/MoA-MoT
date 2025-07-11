import math

# Define the physical parameters for the rope and the lift
# You can change these values to see how the required force changes
l = 10.0  # length of the rope in meters
m = 5.0   # mass of the rope in kilograms
v = 4.0   # speed of the top end of the rope in m/s
g = 9.8   # acceleration due to gravity in m/s^2

print(f"Calculating the required force F for a rope with:")
print(f"  Length (l) = {l} m")
print(f"  Mass (m) = {m} kg")
print(f"  Final Speed (v) = {v} m/s")
print("-" * 30)

# The derived formula for the force F is:
# F = (m * v^2) / (6 * l) + (m * g) / 2

# Calculate the two components of the force
kinetic_energy_term = (m * v**2) / (6 * l)
potential_energy_term = (m * g) / 2

# Calculate the total force
total_force = kinetic_energy_term + potential_energy_term

# Print the final equation with the numbers substituted in, as requested.
print("Final Formula:")
print("F = (m * v^2) / (6 * l) + (m * g) / 2\n")

print("Calculation Steps:")
print(f"F = ({m} * {v}**2) / (6 * {l}) + ({m} * {g}) / 2")
print(f"F = ({m * v**2}) / ({6 * l}) + ({m * g}) / 2")
print(f"F = {kinetic_energy_term:.4f} + {potential_energy_term:.4f}")

# Print the final result
print("-" * 30)
print(f"The exact force F required is: {total_force:.4f} Newtons")