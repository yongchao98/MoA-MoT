# Engine speed in Revolutions Per Minute
engine_rpm = 4500

# ECU processor speed in Megahertz
ecu_mhz = 450

# 1. Convert engine speed from RPM to Revolutions Per Second (RPS)
revolutions_per_second = engine_rpm / 60

# 2. Convert ECU speed from MHz to cycles per second (Hz)
ecu_hz = ecu_mhz * 1_000_000

# 3. Calculate the total number of CPU cycles available per revolution.
# This is the theoretical maximum for the algorithm since the interrupt occurs once per revolution.
cycles_per_revolution = ecu_hz / revolutions_per_second

# Print the final equation and the result
print(f"The theoretical maximum number of cycles is calculated as:")
print(f"({int(ecu_hz)} cycles/sec) / ({int(revolutions_per_second)} revolutions/sec) = {int(cycles_per_revolution)} cycles/revolution")
print(f"Therefore, the ECU can run a maximum of {int(cycles_per_revolution):,} cycles for the algorithm.")

<<<6000000>>>