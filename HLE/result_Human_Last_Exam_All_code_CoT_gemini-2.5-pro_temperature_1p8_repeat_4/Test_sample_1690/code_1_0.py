import math

# Step 1: Define initial values and Wuxing system parameters
P = 100.0
x = 4.0
y = 5.0
z = 1.0
T_half = 400.0

# --- Calculation for a: Camera shutdown time ---
P_req_camera = x + y + z
t_camera = T_half * math.log(P / P_req_camera) / math.log(2)
t_camera_rounded = round(t_camera)

print("a. Time to turn off camera subsystem")
# Using the less optimized formula for clarity as requested by the prompt format
# t = T_half * (ln(P) - ln(P_req)) / ln(2)
print("Equation: t = T_half * (ln(P) - ln(x+y+z)) / ln(2)")
print("Values: t =", int(T_half), "* (ln(", int(P), ") - ln(", int(P_req_camera), ")) / ln(2) =", t_camera_rounded)
print("")

# --- Calculation for b: Sensor shutdown time ---
P_req_sensor = x + z
t_sensor = T_half * math.log(P / P_req_sensor) / math.log(2)
t_sensor_rounded = round(t_sensor)

print("b. Time to turn off sensor subsystem")
print("Equation: t = T_half * (ln(P) - ln(x+z)) / ln(2)")
print("Values: t =", int(T_half), "* (ln(", int(P), ") - ln(", int(P_req_sensor), ")) / ln(2) =", t_sensor_rounded)
print("")

# --- Calculation for c: Memory usage ---
# An efficient program declares P,x,y,z as ints (5D each)
# and t_half, ln2_val, t_camera, t_sensor as fracs (6D each).
mem_int_count = 4
mem_int_size = 5
mem_frac_count = 4 # t_half, ln2, result_a, result_b
mem_frac_size = 6
mem_usage = (mem_int_count * mem_int_size) + (mem_frac_count * mem_frac_size)

print("c. Memory usage in D for program variables")
print("Equation: (int_vars * 5D) + (frac_vars * 6D)")
print("Values: (", mem_int_count, "*", mem_int_size, ") + (", mem_frac_count, "*", mem_frac_size, ") =", mem_usage)
print("")

# --- Calculation for d: Number of exp function calls ---
# The efficient formula requires ln(2), ln(10), and ln(20).
# Number of exp calls for ln(v) is ceil(log2(v) + 9*log2(10))
log2_10_term = 9 * math.log2(10)

val_for_ln1 = 2
val_for_ln2 = P / P_req_camera
val_for_ln3 = P / P_req_sensor

calls_for_ln2 = math.ceil(math.log2(val_for_ln1) + log2_10_term)
calls_for_ln10 = math.ceil(math.log2(val_for_ln2) + log2_10_term)
calls_for_ln20 = math.ceil(math.log2(val_for_ln3) + log2_10_term)
total_exp_calls = calls_for_ln2 + calls_for_ln10 + calls_for_ln20

print("d. Number of 'exp' function calls")
print("Equation: calls for ln(2) + calls for ln(10) + calls for ln(20)")
print("Values:", calls_for_ln2, "+", calls_for_ln10, "+", calls_for_ln20, "=", total_exp_calls)
print("")

# --- Final Answer ---
final_answer = f"{t_camera_rounded}:{t_sensor_rounded}:{mem_usage}:{total_exp_calls}"
print("Final Answer (a:b:c:d):")
print(f"<<<{final_answer}>>>")