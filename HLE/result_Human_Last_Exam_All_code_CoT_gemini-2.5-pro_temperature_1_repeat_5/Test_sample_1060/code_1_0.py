# Define the given parameters
engine_rpm = 4500  # Revolutions Per Minute
ecu_speed_mhz = 450   # Megahertz

# 1. Convert engine speed from RPM to Revolutions Per Second (RPS)
revolutions_per_second = engine_rpm / 60

# 2. Convert ECU speed from MHz to Hz (cycles per second)
ecu_speed_hz = ecu_speed_mhz * 1_000_000

# 3. Calculate the theoretical maximum number of cycles per revolution
# This is the total cycles per second divided by the revolutions per second.
max_cycles_per_revolution = ecu_speed_hz / revolutions_per_second

# --- Output the results ---
print("This script calculates the maximum CPU cycles available for an ECU algorithm per engine revolution.")
print("-" * 80)
print(f"Given Engine Speed: {engine_rpm} RPM")
print(f"Given ECU Processor Speed: {ecu_speed_mhz} MHz")
print("-" * 80)
print("Step 1: Calculate Revolutions Per Second (RPS)")
print(f"{engine_rpm} RPM / 60 seconds/minute = {revolutions_per_second} RPS")
print("\nStep 2: Calculate ECU Speed in Cycles Per Second (Hz)")
print(f"{ecu_speed_mhz} MHz * 1,000,000 Hz/MHz = {int(ecu_speed_hz)} Hz")
print("\nStep 3: Calculate Max Cycles per Revolution (Cycles/RPS)")
print("The final equation is: ECU Speed (Hz) / Engine Speed (RPS)")
print(f"{int(ecu_speed_hz)} / {revolutions_per_second} = {int(max_cycles_per_revolution)}")
print("-" * 80)
print(f"The theoretical maximum number of cycles is: {int(max_cycles_per_revolution)}")
