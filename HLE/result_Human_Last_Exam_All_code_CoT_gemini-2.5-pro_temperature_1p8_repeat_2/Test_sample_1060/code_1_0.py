import math

# --- Inputs ---
# Engine speed in Revolutions Per Minute
rpm = 4500
# ECU processor speed in MegaHertz
cpu_freq_mhz = 450
# Interrupts per revolution
interrupts_per_rev = 1

# --- Calculations ---
# 1. Convert CPU frequency from MHz to Hz
cpu_freq_hz = cpu_freq_mhz * 1_000_000

# 2. Calculate engine speed in Revolutions Per Second (RPS)
revolutions_per_second = rpm / 60

# 3. Calculate the total number of CPU cycles available per revolution
# This is the total cycles per second divided by the revolutions per second.
max_cycles_per_revolution = cpu_freq_hz / revolutions_per_second

# --- Output ---
print(f"To find the theoretical maximum cycles for the algorithm, we follow these steps:")
print(f"1. Determine the number of revolutions per second (RPS):")
print(f"   {rpm} RPM / 60 = {revolutions_per_second} RPS")
print("\n2. Determine the number of CPU cycles per second (Hz):")
print(f"   {cpu_freq_mhz} MHz * 1,000,000 = {int(cpu_freq_hz)} Hz")
print("\n3. Divide cycles per second by revolutions per second to get cycles per revolution:")
print("\nFinal Equation:")
print(f"{int(cpu_freq_hz)} / {revolutions_per_second} = {int(max_cycles_per_revolution)}")
<<<6000000>>>