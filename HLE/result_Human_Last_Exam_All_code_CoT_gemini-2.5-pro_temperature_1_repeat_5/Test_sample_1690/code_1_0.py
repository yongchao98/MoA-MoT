import math

# --- Given Parameters ---
# Initial power P (W)
P = 100
# Power requirements for subsystems (W)
x = 4  # control-and-communication
y = 5  # camera
z = 1  # sensor
# Half-life of the battery (days)
T_half = 400

# --- Part a: Camera subsystem shutdown time ---
# The camera turns off when the available power is no longer sufficient for all three subsystems.
# Power required for all three subsystems:
P_req1 = x + y + z
# Using the derived formula: t = T_half * log2(P_initial / P_required)
t1 = T_half * math.log2(P / P_req1)
# Round to the nearest total day
a = round(t1)

# --- Part b: Sensor subsystem shutdown time ---
# The sensor turns off when the available power is no longer sufficient for the remaining two subsystems.
# Power required for control and sensor subsystems:
P_req2 = x + z
# Using the same formula:
t2 = T_half * math.log2(P / P_req2)
# Round to the nearest total day
b = round(t2)

# --- Part c: Memory usage in D for variables ---
# To be memory-efficient, the C program would reuse variables.
# The necessary variables of type 'frac' are:
# 5 for inputs (P, x, y, z, T_half)
# 1 for a pre-calculated constant (ln2)
# 4 for calculation temporaries (P_req, ratio_inv, ln_ratio_inv, t)
num_frac_vars = 5 + 1 + 4
# The size of one 'frac' variable is 6D (2D+2D+2D).
size_of_frac_in_D = 6
c = num_frac_vars * size_of_frac_in_D

# --- Part d: Number of calls to the 'exp' function ---
# The calculation requires three calls to the ln function: ln(2), ln(10), and ln(20).
# The number of 'exp' calls within ln(val) is ceil(log2(val * 1e9)).
precision_val = 1e9
val1 = 2
val2 = P / P_req1  # This is 10
val3 = P / P_req2  # This is 20

# Calculate the number of calls for each ln invocation
calls_for_ln2 = math.ceil(math.log2(val1 * precision_val))
calls_for_ln10 = math.ceil(math.log2(val2 * precision_val))
calls_for_ln20 = math.ceil(math.log2(val3 * precision_val))

# Total number of exp calls is the sum
d = calls_for_ln2 + calls_for_ln10 + calls_for_ln20

# --- Final Output ---
# Print the answers for a, b, c, and d in the specified format.
print(f"{a}:{b}:{c}:{d}")