import math

# Step 9: We have the equation for p_f: p_f^4 - 18*p_f - 27 = 0.
# Step 10: We find by inspection that p_f = 3 is a solution.
p_f = 3

# Verify that p_f = 3 is a root of the equation.
eq_val = p_f**4 - 18*p_f - 27
print(f"Verifying the root p_f = {p_f}:")
print(f"{p_f}^4 - 18*{p_f} - 27 = {p_f**4} - {18*p_f} - 27 = {eq_val}")
if eq_val == 0:
    print("p_f = 3 is a correct root.\n")
else:
    print("p_f = 3 is not a correct root.\n")

# Step 11 & 12: Calculate x0 using the formula x(p) = -3(p+1)/sqrt(2p+3).
p = p_f
# The denominator is 2*p + 3
denominator_val = 2 * p + 3
# The numerator is -3 * (p + 1)
numerator_val = -3 * (p + 1)

# Calculate x0
x0 = numerator_val / math.sqrt(denominator_val)

print("The final position x0 is calculated using the equation:")
print("x0 = -3 * (p + 1) / sqrt(2*p + 3)")
print("\nSubstituting p = 3:")
print(f"x0 = -3 * ({p} + 1) / sqrt(2*{p} + 3)")
print(f"x0 = {numerator_val} / sqrt({denominator_val})")
print(f"x0 = {numerator_val} / {math.sqrt(denominator_val)}")
print(f"x0 = {x0}")
