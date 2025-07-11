import math

# Step 1: Define initial parameters
P = 100  # Initial power in Watts
x = 4    # Control subsystem power in Watts
y = 5    # Camera subsystem power in Watts
z = 1    # Sensor subsystem power in Watts
T_half = 400 # Half-life in days

# --- Part a: Camera subsystem turn-off time ---
print("--- Part a: Camera turn-off time (t1) ---")
# Equation: P * (1/2)^(t1 / 400) = x + y + z
# Solved for t1: t1 = 400 * ln(P / (x + y + z)) / ln(2)
pwr_threshold_a = x + y + z
ratio_a = P / pwr_threshold_a
ln_ratio_a = math.log(ratio_a)
ln_2 = math.log(2)
t1 = T_half * ln_ratio_a / ln_2
a = round(t1)

print(f"The equation is: t1 = {T_half} * ln({P} / ({x} + {y} + {z})) / ln(2)")
print(f"Intermediate calculation: t1 = {T_half} * ln({ratio_a:.2f}) / ln(2)")
print(f"Using logarithm values: t1 = {T_half} * {ln_ratio_a:.6f} / {ln_2:.6f}")
print(f"Result: t1 = {t1:.2f} days")
print(f"Rounded value (a): {a}\n")


# --- Part b: Sensor subsystem turn-off time ---
print("--- Part b: Sensor turn-off time (t2) ---")
# Equation: P * (1/2)^(t2 / 400) = x + z
# Solved for t2: t2 = 400 * ln(P / (x + z)) / ln(2)
pwr_threshold_b = x + z
ratio_b = P / pwr_threshold_b
ln_ratio_b = math.log(ratio_b)
t2 = T_half * ln_ratio_b / ln_2
b = round(t2)

print(f"The equation is: t2 = {T_half} * ln({P} / ({x} + {z})) / ln(2)")
print(f"Intermediate calculation: t2 = {T_half} * ln({ratio_b:.2f}) / ln(2)")
print(f"Using logarithm values: t2 = {T_half} * {ln_ratio_b:.6f} / {ln_2:.6f}")
print(f"Result: t2 = {t2:.2f} days")
print(f"Rounded value (b): {b}\n")


# --- Part c: Memory usage for variables ---
print("--- Part c: Memory usage in D ---")
# Most efficient program would use:
# 4 int variables (P, x, y, z): 5D each
# 1 frac variable to store ln(2)
# 1 frac variable to store the result 't' (reused for t1 and t2)
# A 'frac' is 3 chars, and each char is 2D. So a frac is 6D.
int_size_D = 5
frac_size_D = 3 * 2 # 3 chars * 2D/char
num_ints = 4
num_fracs = 2
c = (num_ints * int_size_D) + (num_fracs * frac_size_D)

print(f"Calculation: ({num_ints} ints * {int_size_D} D) + ({num_fracs} fracs * {frac_size_D} D)")
print(f"Memory usage (c): {c} D\n")


# --- Part d: Number of calls to exp() ---
print("--- Part d: Total calls to exp() function ---")
# exp() is called inside ln(v) once per iteration.
# Iterations for ln(v) ~ ceil(log2(v / precision)), where precision is 1e-9.
# The program calculates ln(2), ln(10), and ln(20).
precision = 1e-9
def count_exp_calls(value):
  # Number of iterations is ceil(log2(value / precision))
  return math.ceil(math.log2(value / precision))

calls_for_ln2 = count_exp_calls(2)
calls_for_ln10 = count_exp_calls(ratio_a)
calls_for_ln20 = count_exp_calls(ratio_b)
d = calls_for_ln2 + calls_for_ln10 + calls_for_ln20

print("Calls from ln(2): ceil(log2(2 / 1e-9)) = " + str(calls_for_ln2))
print("Calls from ln(10): ceil(log2(10 / 1e-9)) = " + str(calls_for_ln10))
print("Calls from ln(20): ceil(log2(20 / 1e-9)) = " + str(calls_for_ln20))
print(f"Total exp() calls (d): {calls_for_ln2} + {calls_for_ln10} + {calls_for_ln20} = {d}\n")


# --- Final Answer ---
print(f"Final Answer (a:b:c:d): {a}:{b}:{c}:{d}")