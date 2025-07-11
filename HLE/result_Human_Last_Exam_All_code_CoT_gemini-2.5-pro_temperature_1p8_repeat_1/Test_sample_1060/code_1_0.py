# --- Given Parameters ---
# Engine speed in Revolutions Per Minute (RPM)
engine_rpm = 4500

# ECU processor speed in MegaHertz (MHz)
ecu_speed_mhz = 450

# --- Calculations ---

# 1. Convert engine speed from RPM to Revolutions Per Second (RPS)
revolutions_per_second = engine_rpm / 60

# 2. Convert ECU speed from MHz to cycles per second (Hz)
# 1 MHz = 1,000,000 Hz
cycles_per_second = ecu_speed_mhz * 1_000_000

# 3. Calculate the maximum number of cycles available per revolution.
# This is found by dividing the number of cycles per second by the revolutions per second.
max_cycles_per_revolution = cycles_per_second / revolutions_per_second

# --- Output the results ---
print("Calculating the maximum number of ECU cycles per engine revolution.")
print("-" * 60)
print(f"Engine speed: {engine_rpm} RPM")
print(f"ECU processor speed: {ecu_speed_mhz} MHz\n")

print("Step 1: Convert units to a 'per second' basis.")
print(f"ECU Cycles per Second = {ecu_speed_mhz} MHz * 1,000,000 = {int(cycles_per_second)} cycles/sec")
print(f"Engine Revolutions per Second = {engine_rpm} RPM / 60 = {revolutions_per_second} revs/sec\n")

print("Step 2: Calculate the number of cycles per revolution.")
print("This is the final equation:")
print(f"{int(cycles_per_second)} (cycles/sec) / {revolutions_per_second} (revs/sec) = {int(max_cycles_per_revolution)} (cycles/rev)\n")

print(f"The theoretical maximum number of cycles you can run between interrupts is {int(max_cycles_per_revolution)}.")