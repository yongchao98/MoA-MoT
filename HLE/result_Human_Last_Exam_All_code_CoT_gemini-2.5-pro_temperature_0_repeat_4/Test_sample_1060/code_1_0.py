# Define the initial parameters from the problem description
engine_rpm = 4500
ecu_speed_mhz = 450

# Convert engine speed from Revolutions Per Minute (RPM) to Revolutions Per Second (RPS)
# There are 60 seconds in a minute.
engine_rps = engine_rpm / 60

# Convert ECU speed from MegaHertz (MHz) to Hertz (Hz), which is cycles per second
# 1 MHz = 1,000,000 Hz
ecu_speed_hz = ecu_speed_mhz * 1_000_000

# Calculate the maximum number of cycles per revolution by dividing
# the total cycles per second (Hz) by the revolutions per second (RPS).
max_cycles = ecu_speed_hz / engine_rps

print("To find the maximum number of cycles per revolution, we use the following formula:")
print("Max Cycles = (ECU Speed in Hz) / (Engine Speed in RPS)\n")

print("The final calculation is:")
# Print each number in the final equation as requested
print(f"{int(ecu_speed_hz)} / {int(engine_rps)} = {int(max_cycles)}")