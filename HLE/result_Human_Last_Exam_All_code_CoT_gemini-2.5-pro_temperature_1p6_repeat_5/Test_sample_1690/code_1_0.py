import math

# --- Inputs ---
P = 100
x = 4
y = 5
z = 1

print("Solving for P=100, x=4, y=5, z=1\n")

# --- Part a: Camera subsystem turn-off time ---
p_threshold_cam = x + y + z
ratio_cam = P / p_threshold_cam
t_cam = 400 * math.log(ratio_cam) / math.log(2)
t_cam_rounded = round(t_cam)

print("--- Part a: Camera Turn-off Time ---")
print(f"Equation: t = 400 * ln(P / (x + y + z)) / ln(2)")
print(f"Calculation: t = 400 * ln({P} / ({x} + {y} + {z})) / ln(2)")
print(f"Calculation: t = 400 * ln({P} / {p_threshold_cam}) / ln(2)")
print(f"Result: t = {t_cam:.2f} days")
print(f"a = {t_cam_rounded} (rounded to nearest day)\n")


# --- Part b: Sensor subsystem turn-off time ---
p_threshold_sens = x + z
ratio_sens = P / p_threshold_sens
t_sens = 400 * math.log(ratio_sens) / math.log(2)
t_sens_rounded = round(t_sens)

print("--- Part b: Sensor Turn-off Time ---")
print(f"Equation: t = 400 * ln(P / (x + z)) / ln(2)")
print(f"Calculation: t = 400 * ln({P} / ({x} + {z})) / ln(2)")
print(f"Calculation: t = 400 * ln({P} / {p_threshold_sens}) / ln(2)")
print(f"Result: t = {t_sens:.2f} days")
print(f"b = {t_sens_rounded} (rounded to nearest day)\n")


# --- Part c: Memory Usage in D ---
int_d = 5
frac_d = 6 # (signed char (2D) + unsigned char (2D) + signed char (2D))
num_ints = 4 # P, x, y, z
num_fracs = 3 # ln2, ln_ratio, time_days
mem_usage = (num_ints * int_d) + (num_fracs * frac_d)

print("--- Part c: Memory Usage ---")
print(f"Optimal program uses {num_ints} int variables and {num_fracs} frac variables.")
print(f"Equation: Memory = (num_ints * size_of_int) + (num_fracs * size_of_frac)")
print(f"Calculation: Memory = ({num_ints} * {int_d}D) + ({num_fracs} * {frac_d}D)")
print(f"Result: Memory = {num_ints*int_d}D + {num_fracs*frac_d}D")
print(f"c = {mem_usage}D\n")


# --- Part d: Number of exp() calls ---
def get_exp_calls(val):
    # k = ceil(log2(x * 1e9))
    return math.ceil(math.log2(val * 1e9))

calls_ln2 = get_exp_calls(2)
calls_ln10 = get_exp_calls(ratio_cam)
calls_ln20 = get_exp_calls(ratio_sens)
total_calls = calls_ln2 + calls_ln10 + calls_ln20

print("--- Part d: Number of exp() Calls ---")
print("Equation: calls = calls_for_ln(2) + calls_for_ln(ratio_a) + calls_for_ln(ratio_b)")
print(f"Calls for ln(2)  = ceil(log2(2 * 1e9)) = {calls_ln2}")
print(f"Calls for ln({ratio_cam}) = ceil(log2({ratio_cam} * 1e9)) = {calls_ln10}")
print(f"Calls for ln({ratio_sens}) = ceil(log2({ratio_sens} * 1e9)) = {calls_ln20}")
print(f"Total calls = {calls_ln2} + {calls_ln10} + {calls_ln20}")
print(f"d = {total_calls}\n")


# --- Final Answer ---
final_answer = f"{t_cam_rounded}:{t_sens_rounded}:{mem_usage}:{total_calls}"
print(f"Final Answer (a:b:c:d): {final_answer}")