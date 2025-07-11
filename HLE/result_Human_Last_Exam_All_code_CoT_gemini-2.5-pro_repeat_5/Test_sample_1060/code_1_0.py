# Given values
rpm = 4500  # Revolutions Per Minute
cpu_speed_mhz = 450  # ECU processor speed in MegaHertz

# 1. Convert RPM to Revolutions Per Second (RPS)
revolutions_per_second = rpm / 60

# 2. Convert CPU speed from MHz to Hz (cycles per second)
cpu_speed_hz = cpu_speed_mhz * 1_000_000

# 3. Calculate the time it takes for one revolution
# This is the window of time available for the algorithm
# time_per_revolution = 1 / revolutions_per_second

# 4. Calculate the maximum number of cycles available in that time window
# Max Cycles = CPU Speed (cycles/sec) * Time per Revolution (sec/rev)
# This can be simplified to: CPU Speed (cycles/sec) / Revolutions per Second (rev/sec)
max_cycles = cpu_speed_hz / revolutions_per_second

# Print the final equation with all the numbers
print(f"The theoretical maximum number of cycles per revolution is calculated as follows:")
print(f"({int(cpu_speed_hz):,} cycles/second) / ({int(revolutions_per_second):,} revolutions/second) = {int(max_cycles):,} cycles/revolution")

# The final answer
# print(f"\nFinal Answer: {int(max_cycles)}")
<<<6000000>>>