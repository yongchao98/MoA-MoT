import math

# --- Monk's Challenge Parameters ---
# You can change these values to see how the required force changes.

# Mass of the rope in kilograms (m)
m = 5.0

# Length of the rope in meters (l)
l = 10.0

# Final speed of the rope's end in meters per second (v)
v = 2.0

# Acceleration due to gravity in m/s^2 (g)
g = 9.8

# --- Calculation ---
# The formula is derived from the work-energy theorem: F * l = (m*g*l/2) + (1/2*m*v^2)
# Solving for F gives: F = m/2 * (g + v^2/l)

force = (m / 2.0) * (g + (v**2 / l))

# --- Output the Solution ---
print("To find the mystical force F, we use the work-energy theorem.")
print("Work done by force F = Change in Potential Energy + Change in Kinetic Energy")
print(f"W_F = ΔU + ΔK")
print(f"F * l = (m * g * l / 2) + (1/2 * m * v^2)")
print("\nSolving for F, the formula is: F = (m/2) * (g + v^2/l)")

print("\n--- Substituting the given values ---")
print(f"Mass (m) = {m} kg")
print(f"Length (l) = {l} m")
print(f"Speed (v) = {v} m/s")
print(f"Gravity (g) = {g} m/s^2")

print("\n--- The Final Equation ---")
# Printing each number in the final equation as requested
print(f"F = ({m} / 2) * ({g} + {v}**2 / {l})")
print(f"F = {m/2} * ({g} + {v**2 / l})")
print(f"F = {m/2} * {g + (v**2 / l)}")
print(f"\nThe calculated force F is: {force:.2f} Newtons")
