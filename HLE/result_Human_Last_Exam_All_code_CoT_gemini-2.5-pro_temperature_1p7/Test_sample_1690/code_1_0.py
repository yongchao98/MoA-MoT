import math

# --- Input Parameters from the problem description ---
P_initial = 100.0   # Initial power in Watts
x = 4.0             # Control subsystem power
y = 5.0             # Camera subsystem power
z = 1.0             # Sensor subsystem power
T_half = 400.0      # Power half-life in days

# --- a. Camera subsystem shutdown time ---
# Shutdown occurs when power drops below the requirement for all systems.
power_threshold_cam = x + y + z
# The governing equation is t = T_half * log2(P_initial / P_target).
t_cam = T_half * math.log2(P_initial / power_threshold_cam)
answer_a = round(t_cam)
print(f"a. The camera turns off when power is below {power_threshold_cam}W. The equation is t = {T_half} * log2({P_initial}/{power_threshold_cam}). This happens after {answer_a} days.")


# --- b. Sensor subsystem shutdown time ---
# Shutdown occurs when power drops below the requirement for control and sensor systems.
power_threshold_sen = x + z
# Since power_threshold_sen (5W) is half of power_threshold_cam (10W),
# this will occur one half-life after the camera shutdown.
t_sen = t_cam + T_half
answer_b = round(t_sen)
print(f"b. The sensor turns off when power is below {power_threshold_sen}W. The equation is t = {answer_a} + {T_half}. This happens after {answer_b} days.")


# --- c. Memory usage in D ---
# An efficient Wuxing C program requires a minimum of 3 'frac' variables.
# Each 'frac' variable has a size of (2D + 2D + 2D) = 6D.
frac_size_D = 6
min_frac_variables = 3
answer_c = frac_size_D * min_frac_variables
print(f"c. The most memory-efficient program uses {min_frac_variables} 'frac' variables, each 6D in size. The equation is memory = {min_frac_variables} * {frac_size_D}. Total usage is {answer_c} D.")


# --- d. Number of 'exp' function calls ---
# The efficient program calculates ln(2) and ln(10).
# Number of calls for ln(v) is ceil(log2(v * 10^9)).
# Number of calls for ln(2)
calls_ln2 = math.ceil(1 + 9 * math.log2(10))
# Number of calls for ln(10)
calls_ln10 = math.ceil(10 * math.log2(10))
answer_d = calls_ln2 + calls_ln10
print(f"d. The program calls ln(2) and ln(10). The equation is total_calls = calls_for_ln(2) + calls_for_ln(10) = {calls_ln2} + {calls_ln10}. Total expensive calls are {answer_d}.")


# --- Final Answer in a:b:c:d format ---
final_answer = f"{answer_a}:{answer_b}:{answer_c}:{answer_d}"
print(f"\nFinal Answer: {final_answer}")