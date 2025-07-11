import math

# Define input parameters based on the problem description
P_initial = 100.0
x = 4.0
y = 5.0
z = 1.0
T_half = 400.0

# a. When it turn off the camera subsystem (rounded to total days)?
# This occurs when power drops below the requirement for all three systems.
P_threshold_a = x + y + z
# The formula for time is t = T_half * ln(P_initial / P_threshold) / ln(2)
Ratio_a = P_initial / P_threshold_a
t_a = T_half * math.log(Ratio_a) / math.log(2)
a = round(t_a)

# b. When it turn off the sensor subsystem (rounded to total days)?
# This occurs when power drops below the requirement for the control and sensor systems.
P_threshold_b = x + z
Ratio_b = P_initial / P_threshold_b
t_b = T_half * math.log(Ratio_b) / math.log(2)
b = round(t_b)

# c. What is the memory usage in D for variables of this program?
# An efficient C program would have 7 main `frac` variables (P,x,y,z,T_half,ln2,result).
# frac = 3 * char = 3 * 2D = 6D. Memory_main = 7 * 6D = 42D.
# The `ln` function uses 4 local `frac` variables: 4 * 6D = 24D.
# The `exp` function uses 2 local `frac`s and 1 `char`: 2 * 6D + 2D = 14D.
# Peak memory usage is the sum during the deepest call (main->ln->exp).
mem_main = 7 * 6
mem_ln = 4 * 6
mem_exp = 2 * 6 + 2
c = mem_main + mem_ln + mem_exp

# d. What is the number of time this program call function exp?
# The program calls ln(2), ln(Ratio_a=10), and ln(Ratio_b=20).
# The number of exp calls for ln(N) is ceil(log2(N * 10^9)).
def calculate_exp_calls(N):
    return math.ceil(math.log2(N * 1e9))

calls_for_ln2 = calculate_exp_calls(2)
calls_for_ln10 = calculate_exp_calls(Ratio_a)
calls_for_ln20 = calculate_exp_calls(Ratio_b)
d = calls_for_ln2 + calls_for_ln10 + calls_for_ln20

# Print the final result in the required format a:b:c:d
# The prompt also asks to show the numbers in the final equation.
print("Calculation for 'a' (Camera Shutdown):")
print(f"t = {T_half} * ln({P_initial} / ({x} + {y} + {z})) / ln(2)")
print(f"t = {T_half} * ln({Ratio_a}) / ln(2) = {t_a:.2f} days")
print(f"Rounded days = {a}")
print("\nCalculation for 'b' (Sensor Shutdown):")
print(f"t = {T_half} * ln({P_initial} / ({x} + {z})) / ln(2)")
print(f"t = {T_half} * ln({Ratio_b}) / ln(2) = {t_b:.2f} days")
print(f"Rounded days = {b}")
print(f"\nCalculation for 'c' (Memory Usage): {c}D")
print(f"Calculation for 'd' (exp calls): {d}")

print(f"\n<<<1329:1729:80:100>>>")