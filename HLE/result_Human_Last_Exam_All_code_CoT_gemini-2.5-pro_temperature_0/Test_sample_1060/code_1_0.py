# Given values
engine_rpm = 4500  # Revolutions Per Minute
ecu_mhz = 450      # MegaHertz

# 1. Convert engine speed to Revolutions Per Second (RPS)
# There are 60 seconds in a minute.
engine_rps = engine_rpm / 60

# 2. Convert ECU speed to cycles per second (Hz)
# 1 MHz = 1,000,000 Hz
ecu_hz = ecu_mhz * 1_000_000

# 3. Calculate the maximum number of cycles per revolution (per interrupt)
# This is the total cycles per second divided by the revolutions per second.
max_cycles_per_revolution = ecu_hz / engine_rps

# 4. Print the final equation and the result
print(f"The engine runs at {engine_rpm} RPM, which is {engine_rps:.0f} revolutions per second.")
print(f"The ECU runs at {ecu_mhz} MHz, which is {ecu_hz:,.0f} cycles per second.")
print("\nTo find the maximum cycles per revolution, we divide the ECU cycles per second by the engine revolutions per second.")
print("Final Equation:")
print(f"{ecu_hz:,.0f} cycles/sec / {engine_rps:.0f} revs/sec = {max_cycles_per_revolution:,.0f} cycles/rev")

print(f"\nThe theoretical maximum number of cycles is {int(max_cycles_per_revolution):,}.")