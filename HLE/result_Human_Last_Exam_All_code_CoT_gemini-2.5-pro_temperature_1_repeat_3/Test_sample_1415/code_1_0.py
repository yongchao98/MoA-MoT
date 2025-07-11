import math

# The integral of F3'(x)/F3(x) from 0 to 1 is ln(F3(1)) - ln(F3(0)).

# Calculate F3(1)
# F3(1) = ln(1 + e^(sin(1)))
sin_1 = math.sin(1)
exp_sin_1 = math.exp(sin_1)
f3_at_1 = math.log(1 + exp_sin_1)

# Calculate F3(0)
# F3(0) = ln(1 + e^(sin(0))) = ln(1 + e^0) = ln(2)
f3_at_0 = math.log(2)

print(f"The final equation for the integral V is: ln(F3(1)) - ln(F3(0))")
print(f"The value of F3(1) is: {f3_at_1}")
print(f"The value of F3(0) is: {f3_at_0}")


# Calculate the value of the integral V
# V = ln(F3(1)) - ln(F3(0))
# Note: ln(a) - ln(b) = ln(a/b)
V = math.log(f3_at_1) - math.log(f3_at_0)

# Calculate the final result as requested
final_value = 10000 * V
closest_integer = round(final_value)

print(f"The value of the integral V is: {V}")
print(f"The value of 10000 * V is: {final_value}")
print(f"The closest integer to 10000 * V is: {closest_integer}")

# The final answer in the required format
# print(f"\n<<<{closest_integer}>>>")