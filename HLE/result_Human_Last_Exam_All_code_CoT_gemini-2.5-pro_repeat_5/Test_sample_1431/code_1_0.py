import math

# --- User-defined variables for the rope ---
# You can change these values to match a specific problem.
m = 10.0  # mass of the rope in kg
l = 20.0  # length of the rope in meters
v = 5.0   # final speed of the rope in m/s

# --- Physical constant ---
g = 9.8   # acceleration due to gravity in m/s^2

# --- Explanation of the formula ---
print("To find the force F, we use the Work-Energy Theorem: W_net = ΔK")
print("Work done by applied force (F * l) + Work done by gravity (-m*g*l/2) = Change in Kinetic Energy (0.5*m*v^2)")
print("This gives the formula: F = (m * g / 2) + (m * v^2 / (2 * l))\n")

# --- Calculation ---
# Calculate the component of the force needed to overcome gravity
force_gravity_component = (m * g) / 2
# Calculate the component of the force needed to provide kinetic energy
force_kinetic_component = (m * v**2) / (2 * l)
# The total force is the sum of the two components
F_total = force_gravity_component + force_kinetic_component

# --- Output the result ---
print("--- Calculating the Force F ---")
print(f"Given values:")
print(f"  Mass (m) = {m} kg")
print(f"  Length (l) = {l} m")
print(f"  Final Speed (v) = {v} m/s")
print(f"  Gravity (g) = {g} m/s^2\n")

print("Plugging the values into the formula:")
print(f"F = ({m} * {g} / 2) + ({m} * {v}**2 / (2 * {l}))")
print(f"F = {force_gravity_component} + {force_kinetic_component}")
print(f"\nThe total force F required is: {F_total} Newtons")

print(f"\n<<<{F_total}>>>")