# Given parameters
engine_rpm = 4500  # Revolutions Per Minute
ecu_mhz = 450     # MegaHertz

# --- Step 1: Convert Engine Speed from RPM to RPS ---
# There are 60 seconds in a minute.
engine_rps = engine_rpm / 60
print(f"Engine speed: {engine_rpm} RPM")
print(f"Engine speed in RPS: {engine_rpm} / 60 = {engine_rps:.2f} Revolutions Per Second\n")

# --- Step 2: Calculate the time for one revolution ---
# This is the inverse of the RPS.
time_per_revolution = 1 / engine_rps
print(f"Time per revolution: 1 / {engine_rps:.2f} = {time_per_revolution:.6f} seconds\n")

# --- Step 3: Convert ECU speed from MHz to Cycles Per Second ---
# 1 MHz = 1,000,000 Hz (cycles per second)
ecu_hz = ecu_mhz * 1_000_000
print(f"ECU speed: {ecu_mhz} MHz")
print(f"ECU speed in Hz: {ecu_mhz} * 1,000,000 = {int(ecu_hz)} Cycles Per Second\n")

# --- Step 4: Calculate the theoretical maximum number of cycles per revolution ---
# This is the number of cycles per second multiplied by the seconds per revolution.
max_cycles_per_revolution = ecu_hz * time_per_revolution

# --- Final Answer ---
print("Final Calculation:")
print(f"Max cycles = (ECU Speed in Hz) * (Time Per Revolution)")
print(f"Max cycles = {int(ecu_hz)} * {time_per_revolution:.6f}")
print(f"The theoretical maximum number of cycles is: {int(max_cycles_per_revolution)}")