# Given values
engine_rpm = 4500
ecu_speed_mhz = 450

# 1. Convert engine RPM to Revolutions Per Second (RPS)
# There are 60 seconds in a minute.
engine_rps = engine_rpm / 60

# 2. Convert ECU speed from MHz to Hz (cycles per second)
# 1 MHz = 1,000,000 Hz
ecu_speed_hz = ecu_speed_mhz * 1_000_000

# 3. Calculate the maximum number of cycles per revolution.
# This is the processor's speed in cycles/sec divided by the engine's speed in revs/sec.
# The result is the number of cycles the processor can execute in one revolution.
max_cycles_per_revolution = ecu_speed_hz / engine_rps

# Print the explanation and the final equation with the calculated numbers
print(f"To find the maximum cycles per revolution, we divide the ECU's speed in cycles per second by the engine's speed in revolutions per second.")
print("\nFinal Equation:")
# Using integers for cleaner output, as the results are whole numbers.
print(f"Max Cycles = {int(ecu_speed_hz)} / {int(engine_rps)}")
print(f"Result: {int(max_cycles_per_revolution)} cycles")

<<<6000000>>>