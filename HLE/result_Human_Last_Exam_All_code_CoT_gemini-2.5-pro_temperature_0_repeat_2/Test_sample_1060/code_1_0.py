# Define the given parameters
engine_rpm = 4500  # Revolutions Per Minute
ecu_speed_mhz = 450  # ECU speed in MegaHertz

# Step 1: Convert engine speed to Revolutions Per Second (RPS)
# There are 60 seconds in a minute.
engine_rps = engine_rpm / 60

# Step 2: Convert ECU speed to cycles per second (Hz)
# 1 MHz = 1,000,000 Hz
ecu_speed_hz = ecu_speed_mhz * 1_000_000

# Step 3: Calculate the maximum number of cycles per revolution
# This is the total cycles per second divided by the revolutions per second.
max_cycles_per_revolution = ecu_speed_hz / engine_rps

# Print the final equation with all the numbers
print("The theoretical maximum number of cycles is calculated by dividing the ECU's cycles per second by the engine's revolutions per second.")
print("Equation: ECU Speed (Hz) / Engine Speed (RPS) = Max Cycles per Revolution")
print(f"{int(ecu_speed_hz)} / {int(engine_rps)} = {int(max_cycles_per_revolution)}")