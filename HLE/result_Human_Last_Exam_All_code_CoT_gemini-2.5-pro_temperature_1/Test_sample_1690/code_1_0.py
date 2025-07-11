import math

# --- Inputs and Constants ---
P = 100.0  # Initial power
x = 4.0    # Control subsystem power
y = 5.0    # Camera subsystem power
z = 1.0    # Sensor subsystem power
T_HALF = 400 # Power half-life in days

# --- Part a: Camera Shutdown Time ---
# Power level when camera is turned off
power_needed_a = x + y + z
# Ratio of initial power to needed power
ratio_a = P / power_needed_a
# Calculate time t = 400 * ln(ratio) / ln(2)
time_a = T_HALF * math.log(ratio_a) / math.log(2)
# The first number in the equation is the rounded time for camera shutdown
a_camera_shutdown_days = round(time_a)

# --- Part b: Sensor Shutdown Time ---
# Power level when sensor is turned off
power_needed_b = x + z
# Ratio of initial power to needed power
ratio_b = P / power_needed_b
# Calculate time t = 400 * ln(ratio) / ln(2)
time_b = T_HALF * math.log(ratio_b) / math.log(2)
# The second number in the equation is the rounded time for sensor shutdown
b_sensor_shutdown_days = round(time_b)

# --- Part c: Memory Usage ---
# The most efficient program needs 6 frac variables (P, x, y, z, ln2, result).
# Each frac is 3 * char, and each char is 2D.
# So, size of one frac is 3 * 2D = 6D.
num_variables = 6
size_of_frac_in_D = 6
# The third number in the equation is the total memory usage in D
c_memory_usage_D = num_variables * size_of_frac_in_D

# --- Part d: Number of exp() Calls ---
# Number of calls N for ln(val) is ceil(log2(val) + 9 * log2(10))
# Total calls = calls for ln(2) + calls for ln(ratio_a) + calls for ln(ratio_b)
log2_10 = math.log2(10)
# Calls for ln(2)
calls_for_ln2 = math.ceil(math.log2(2) + 9 * log2_10)
# Calls for ln(10)
calls_for_ln10 = math.ceil(math.log2(10) + 9 * log2_10)
# Calls for ln(20)
calls_for_ln20 = math.ceil(math.log2(20) + 9 * log2_10)
# The fourth number in the equation is the total number of exp() calls
d_exp_calls = calls_for_ln2 + calls_for_ln10 + calls_for_ln20

# --- Final Answer ---
# Print the final answer in the format a:b:c:d
print(f"{a_camera_shutdown_days}:{b_sensor_shutdown_days}:{c_memory_usage_D}:{d_exp_calls}")