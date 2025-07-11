# Given values
engine_rpm = 4500  # revolutions per minute
ecu_speed_mhz = 450 # MegaHertz

# 1. Convert ECU speed from MHz to Hz (cycles per second)
ecu_speed_hz = ecu_speed_mhz * 1_000_000

# 2. Convert engine speed from RPM to RPS (revolutions per second)
revolutions_per_second = engine_rpm / 60.0

# 3. Calculate the maximum number of cycles per revolution
# This is done by dividing the total cycles per second by the revolutions per second.
# Cycles per revolution = (Cycles / Second) / (Revolutions / Second) = Cycles / Revolution
max_cycles_per_revolution = ecu_speed_hz / revolutions_per_second

# Print the explanation and the result
print(f"To find the theoretical maximum number of cycles per revolution, we divide the ECU's cycles per second by the engine's revolutions per second.")
print(f"The equation is:")
print(f"({ecu_speed_hz:,} cycles/sec) / ({revolutions_per_second:.2f} revs/sec) = {max_cycles_per_revolution:,.0f} cycles/revolution")

# The final result
print("\nTheoretical maximum number of cycles per revolution:")
print(f"{int(max_cycles_per_revolution):,}")