# Given values
engine_rpm = 4500  # Revolutions Per Minute
ecu_speed_mhz = 450   # MegaHertz

# 1. Convert engine speed from RPM to Revolutions Per Second (RPS)
revolutions_per_second = engine_rpm / 60.0

# 2. Calculate the time available per revolution (which is one interrupt)
time_per_revolution_s = 1.0 / revolutions_per_second

# 3. Convert ECU speed from MHz to Hz (cycles per second)
ecu_speed_hz = ecu_speed_mhz * 1_000_000

# 4. Calculate the maximum number of cycles available between interrupts
max_cycles = ecu_speed_hz * time_per_revolution_s

# Print the final result and the equation
print(f"Engine RPM: {engine_rpm}")
print(f"ECU Speed: {ecu_speed_mhz} MHz")
print("\nTo find the maximum number of cycles available between interrupts, we use the following equation:")
print("Max Cycles = (ECU Speed in Hz) / (Engine Speed in RPS)")
print(f"Max Cycles = ({int(ecu_speed_hz)}) / ({int(revolutions_per_second)})")
print(f"The theoretical maximum number of cycles is: {int(max_cycles)}")
<<<6000000>>>