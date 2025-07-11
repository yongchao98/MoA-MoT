# Given values
engine_rpm = 4500  # Revolutions Per Minute
ecu_speed_mhz = 450   # MegaHertz

# Step 1: Convert engine RPM to Revolutions Per Second (RPS)
# There are 60 seconds in a minute.
revolutions_per_second = engine_rpm / 60

# Step 2: The interrupt frequency is one per revolution, so it equals the RPS.
interrupt_frequency_hz = revolutions_per_second

# Step 3: Convert ECU speed from MHz to Hz
# 1 MHz = 1,000,000 Hz
ecu_speed_hz = ecu_speed_mhz * 1_000_000

# Step 4: Calculate the theoretical maximum number of cycles per interrupt
# This is the total cycles per second divided by the interrupts per second.
max_cycles = ecu_speed_hz / interrupt_frequency_hz

# Output the final result showing the numbers used in the equation
print(f"To find the maximum cycles per interrupt, we divide the ECU's cycles per second by the engine's interrupts per second.")
print(f"Final Calculation: ({int(ecu_speed_hz):,} cycles/sec) / ({int(interrupt_frequency_hz)} interrupts/sec) = {int(max_cycles):,} cycles per interrupt.")
print(f"\nThe theoretical maximum number of cycles you can run is {int(max_cycles):,}.")
