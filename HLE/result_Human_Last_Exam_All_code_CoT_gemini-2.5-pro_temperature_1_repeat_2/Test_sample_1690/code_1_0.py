import math

# Step 1: Define constants and inputs from the problem description.
P = 100.0
x = 4.0
y = 5.0
z = 1.0
T_half = 400.0

# Step 2: Calculate 'a', the time to turn off the camera subsystem.
# This happens when power is no longer sufficient for all three systems.
P_needed_camera = x + y + z
# Using the formula: t = T_half * ln(P_initial / P_final) / ln(2)
t_camera = T_half * math.log(P / P_needed_camera) / math.log(2)
a = int(round(t_camera))

# Step 3: Calculate 'b', the time to turn off the sensor subsystem.
# This happens when power is no longer sufficient for the control and sensor systems.
P_needed_sensor = x + z
t_sensor = T_half * math.log(P / P_needed_sensor) / math.log(2)
b = int(round(t_sensor))

# Step 4: Calculate 'c', the peak memory usage in Decimal digits (D).
# Data type sizes in D: frac=6, int=5, char=2.
# Memory for an optimized main function's variables:
# 5 frac (P,x,y,z,T_half) + 1 frac (ln2) + 3 frac (reusable temps) + 1 int (result)
mem_main = 5 * 6 + 1 * 6 + 3 * 6 + 1 * 5 # = 59 D
# Memory for local variables in the ln() function: 4 frac = 4 * 6 = 24 D
mem_ln_locals = 24
# Memory for local variables in the exp() function: 2 frac + 1 char = 2 * 6 + 2 = 14 D
mem_exp_locals = 14
# Peak memory is the sum during the deepest call stack: main -> ln -> exp
c = mem_main + mem_ln_locals + mem_exp_locals

# Step 5: Calculate 'd', the total number of calls to the expensive 'exp' function.
# The program calculates ln(2), ln(10), and ln(20).
# Each ln(v) call requires ceil(log2(v * 1e9)) calls to exp().
exp_calls_for_ln2 = math.ceil(math.log2(2 * 1e9))
exp_calls_for_ln10 = math.ceil(math.log2(10 * 1e9))
exp_calls_for_ln20 = math.ceil(math.log2(20 * 1e9))
d = exp_calls_for_ln2 + exp_calls_for_ln10 + exp_calls_for_ln20

# Step 6: Print the final answer in the format a:b:c:d.
print(f"{a}:{b}:{c}:{d}")