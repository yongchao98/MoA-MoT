# --- Problem Variables ---
engine_rpm = 4500  # Revolutions Per Minute
ecu_speed_mhz = 450   # Processor speed in MegaHertz

# --- Plan Step 1: Calculate Revolutions Per Second (RPS) ---
# There are 60 seconds in a minute.
revolutions_per_second = engine_rpm / 60.0

# --- Plan Step 2: Calculate Processor Cycles Per Second (Hz) ---
# 1 MHz is 1,000,000 Hz (cycles per second).
ecu_speed_hz = ecu_speed_mhz * 1_000_000

# --- Plan Step 3: Determine Cycles Per Revolution ---
# This is the number of cycles available between interrupts.
# Formula: (Cycles / Second) / (Revolutions / Second) = Cycles / Revolution
max_cycles_per_revolution = ecu_speed_hz / revolutions_per_second

# --- Output the Results ---
print("This script calculates the theoretical maximum CPU cycles available per engine revolution.")
print(f"Given an engine speed of {engine_rpm} RPM and an ECU speed of {ecu_speed_mhz} MHz.\n")

print("The final calculation is (ECU Speed in Hz) / (Engine Speed in RPS).\n")

print("The final equation with the numbers is:")
# The format specifier ':,' adds thousands separators for readability.
print(f"{int(ecu_speed_hz):,} / {int(revolutions_per_second)} = {int(max_cycles_per_revolution):,}\n")

print(f"Therefore, the theoretical maximum number of cycles you can run in the algorithm is {int(max_cycles_per_revolution):,}.")

print(f"\n<<<{int(max_cycles_per_revolution)}>>>")