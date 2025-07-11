import math

# Based on the analysis, the problem reduces to solving the inequality:
# 1 <= 2 * (1 - 10**-M)
# where M is a positive integer.

# We will solve this step-by-step and print the values at each stage.

print("The final inequality to solve for M is derived from W_max <= 2*C, where W_max=1 and C=1-10**-M.")
print("The inequality is:")

ineq_lhs = 1
ineq_rhs_factor = 2
print(f"{ineq_lhs} <= {ineq_rhs_factor} * (1 - 10**-M)")

print("\nStep 1: Divide both sides by 2.")
val1 = ineq_lhs / ineq_rhs_factor
print(f"{val1} <= 1 - 10**-M")

print("\nStep 2: Rearrange the terms to isolate 10**-M.")
# 10**-M <= 1 - val1
val2 = 1 - val1
print(f"10**-M <= {val2}")

print("\nStep 3: Take the reciprocal of both sides, which reverses the inequality sign.")
# 10**M >= 1 / val2
val3 = 1 / val2
print(f"10**M >= {val3}")

print("\nStep 4: Take the base-10 logarithm of both sides.")
log_val = math.log10(val3)
print(f"M >= log10({val3})")
print(f"M >= {log_val}")

# The problem requires M to be the smallest positive integer.
# We need to find the smallest positive integer M such that M >= log_val.
# The value of log_val is approximately 0.301.
# The smallest integer satisfying this is 1.
final_M = math.ceil(log_val)
if final_M < 1:
    final_M = 1

print(f"\nSince M must be a positive integer, the smallest value for M is {final_M}.")
print("<<<1>>>")