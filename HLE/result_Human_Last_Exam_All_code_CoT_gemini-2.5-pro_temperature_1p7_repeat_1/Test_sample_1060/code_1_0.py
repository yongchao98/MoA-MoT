# Define the given parameters
rpm = 4500  # Engine speed in Revolutions Per Minute
ecu_speed_mhz = 450  # ECU processor speed in Megahertz

# 1. Convert RPM to Revolutions Per Second (RPS)
revolutions_per_second = rpm / 60

# 2. Convert ECU speed from MHz to Hz (cycles per second)
ecu_speed_hz = ecu_speed_mhz * 1_000_000

# 3. Calculate the maximum number of cycles per revolution
# This is done by dividing the ECU's cycles per second by the engine's revolutions per second.
# (Cycles / Second) / (Revolutions / Second) = Cycles / Revolution
max_cycles_per_revolution = ecu_speed_hz / revolutions_per_second

# Print the breakdown of the calculation
print(f"Engine speed: {rpm} RPM")
print(f"ECU processor speed: {ecu_speed_mhz} MHz\n")

print("Step 1: Calculate Revolutions Per Second (RPS)")
print(f"{rpm} RPM / 60 s/min = {revolutions_per_second} RPS\n")

print("Step 2: Calculate ECU Speed in Cycles Per Second (Hz)")
print(f"{ecu_speed_mhz} MHz * 1,000,000 = {int(ecu_speed_hz)} Hz\n")

print("Step 3: Calculate Maximum Cycles Per Revolution")
print("This is the number of processor cycles available between interrupts.")
print(f"Final Equation: {int(ecu_speed_hz)} cycles/sec / {revolutions_per_second} revs/sec")
print(f"Result: {int(max_cycles_per_revolution)} cycles\n")

print(f"The theoretical maximum number of cycles you can run in the algorithm is {int(max_cycles_per_revolution)}.")
<<<6000000>>>