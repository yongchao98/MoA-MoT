# Given values
engine_rpm = 4500  # revolutions per minute
processor_speed_mhz = 450  # megahertz

# 1. Convert engine RPM to Revolutions Per Second (RPS)
# There are 60 seconds in a minute.
revolutions_per_second = engine_rpm / 60

# 2. Convert processor speed from MHz to Hz
# 1 MHz = 1,000,000 Hz
processor_speed_hz = processor_speed_mhz * 1_000_000

# 3. Calculate the maximum number of cycles per revolution (interrupt)
# This is the total cycles per second divided by the revolutions per second.
max_cycles_per_revolution = processor_speed_hz / revolutions_per_second

# --- Output the results ---
print("--- Calculation Steps ---")
print(f"Engine Speed: {engine_rpm} RPM")
print(f"ECU Processor Speed: {processor_speed_mhz} MHz")
print("\nStep 1: Calculate Revolutions Per Second (RPS)")
print(f"{engine_rpm} RPM / 60 seconds/minute = {revolutions_per_second:.2f} RPS")

print("\nStep 2: Calculate Processor Speed in Cycles per Second (Hz)")
print(f"{processor_speed_mhz} MHz * 1,000,000 = {int(processor_speed_hz):,} Hz")

print("\nStep 3: Calculate Maximum Cycles per Revolution")
print("This is the number of available CPU cycles between each engine revolution.")
print(f"Final Equation: {int(processor_speed_hz):,} cycles/sec / {revolutions_per_second:.2f} revolutions/sec")

print(f"\nTheoretical Maximum Cycles = {int(max_cycles_per_revolution):,}")
print("<<<6000000>>>")