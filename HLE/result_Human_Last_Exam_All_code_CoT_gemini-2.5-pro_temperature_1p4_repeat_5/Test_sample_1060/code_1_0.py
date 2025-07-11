# Define the given parameters
engine_rpm = 4500  # Revolutions Per Minute
ecu_speed_mhz = 450  # MegaHertz

# Step 1: Convert Engine RPM to Revolutions Per Second (RPS)
seconds_in_a_minute = 60
engine_rps = engine_rpm / seconds_in_a_minute

# Step 2: Convert ECU Speed from MHz to Hertz (cycles per second)
cycles_per_megahertz = 1_000_000
ecu_speed_hz = ecu_speed_mhz * cycles_per_megahertz

# Step 3: Calculate the maximum number of cycles per revolution (per interrupt)
# This is the total cycles per second divided by the revolutions (interrupts) per second.
max_cycles_per_revolution = ecu_speed_hz / engine_rps

# Output the final equation with all the numbers
print(f"The theoretical maximum number of cycles is calculated by dividing the ECU's cycles per second by the engine's revolutions per second.")
print(f"Calculation: (({ecu_speed_mhz} MHz * {cycles_per_megahertz}) cycles/sec) / (({engine_rpm} RPM / {seconds_in_a_minute}) revs/sec)")
print(f"Result: {int(ecu_speed_hz)} cycles/sec / {int(engine_rps)} revs/sec = {int(max_cycles_per_revolution)} cycles/rev")

print("\nTherefore, the theoretical maximum number of cycles available for the algorithm is:")
print(int(max_cycles_per_revolution))
<<<6000000>>>