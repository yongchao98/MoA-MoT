import math

# Plan:
# 1. Calculate time 'a' when the camera is turned off (P_target = x+y+z).
#    The formula is t = half_life * log2(P_initial / P_target).
# 2. Calculate time 'b' when the sensor is turned off (P_target = x+z).
# 3. Calculate memory usage 'c' by summing the variable sizes at maximum
#    call stack depth (main -> ln -> exp).
# 4. Calculate the number of `exp` calls 'd' made by the two necessary `ln` calls
#    (ln(2) and ln(10)).
# 5. Print the final results 'a', 'b', 'c', and 'd' in the specified format.

# --- Parameters based on the problem description ---
P_initial = 100.0
x = 4.0
y = 5.0
z = 1.0
half_life = 400.0

# --- Part a & b: Shutdown Times ---

# a) Camera shutdown time
power_target_a = x + y + z
t_a_float = half_life * math.log2(P_initial / power_target_a)
a = round(t_a_float)

# b) Sensor shutdown time
power_target_b = x + z
t_b_float = half_life * math.log2(P_initial / power_target_b)
b = round(t_b_float)

# --- Part c: Memory Usage ---

# Wuxing data type sizes in decimal digits (D)
size_char = 2
size_int = 5
size_frac = 3 * size_char  # struct frac { char, char, char }

# Memory usage is calculated at maximum stack depth: main() -> ln() -> exp()
# - main() variables: 10 frac, 2 int
# - ln() variables: 5 frac (1 arg + 4 locals)
# - exp() variables: 3 frac (1 arg + 2 locals), 1 char (local)
mem_main = 10 * size_frac + 2 * size_int
mem_ln = 5 * size_frac
mem_exp = 3 * size_frac + 1 * size_char
c = mem_main + mem_ln + mem_exp

# --- Part d: `exp` Function Calls ---

# The most efficient program calculates ln(2) and ln(10).
# The number of `exp` calls is the total iterations in the `ln` function's loop.
# Iterations `k` for `ln(v)` is calculated as `ceil(log2(v * 1e9))`.
calls_for_ln2 = math.ceil(math.log2(2 * 1e9))
calls_for_ln10 = math.ceil(math.log2(10 * 1e9))
d = calls_for_ln2 + calls_for_ln10

# --- Final Output ---
# Print each number for the final answer, formatted as a:b:c:d.
print(f"a: {a}")
print(f"b: {b}")
print(f"c: {c}")
print(f"d: {d}")
print(f"Final Answer: {a}:{b}:{c}:{d}")