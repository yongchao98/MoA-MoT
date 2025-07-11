# Engine speed in Revolutions Per Minute
rpm = 4500

# ECU processor speed in MegaHertz
ecu_mhz = 450

# --- Step 1: Convert RPM to Revolutions Per Second (RPS) ---
# There are 60 seconds in a minute
seconds_per_minute = 60
rps = rpm / seconds_per_minute

# --- Step 2: Convert ECU speed from MHz to Hz ---
# 1 MHz = 1,000,000 Hz
ecu_hz = ecu_mhz * 1_000_000

# --- Step 3: Calculate the maximum cycles per revolution ---
# This is the number of cycles per second divided by the revolutions per second
max_cycles_per_revolution = ecu_hz / rps

# --- Output the final result with the equation ---
# Format numbers with commas for readability
rpm_str = f"{rpm:,}"
ecu_hz_str = f"{ecu_hz:,}"
seconds_str = f"{seconds_per_minute:,}"
result_str = f"{int(max_cycles_per_revolution):,}"

print(f"The calculation is: ({ecu_hz_str} cycles/second) / ({rpm_str} revolutions/minute / {seconds_str} seconds/minute) = {result_str} cycles/revolution")
print(f"The theoretical maximum number of cycles available per revolution is {result_str}.")

# Return the numerical answer in the specified format
# The numerical answer is max_cycles_per_revolution
final_answer = int(max_cycles_per_revolution)