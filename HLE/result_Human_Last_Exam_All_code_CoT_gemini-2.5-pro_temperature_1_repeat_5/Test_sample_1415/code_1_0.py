import math

# Step 1: Define the function F3(x) and calculate its values at x=0 and x=1
# F3(x) = ln(1 + e^(sin(x)))

# Calculate F3(0)
# sin(0) = 0
# e^0 = 1
# 1 + 1 = 2
# ln(2)
f3_0 = math.log(2)

# Calculate F3(1)
# sin(1) is sin(1 radian)
sin_1 = math.sin(1)
exp_sin_1 = math.exp(sin_1)
f3_1 = math.log(1 + exp_sin_1)

# Step 2: Calculate the value of the integral V
# V = ln(F3(1)) - ln(F3(0))
V = math.log(f3_1) - math.log(f3_0)

# Step 3: Calculate the final result
final_value = 10000 * V
closest_integer = round(final_value)

# Print the results as requested
print("The integral evaluates to ln(F₃(1)) - ln(F₃(0))")
print(f"F₃(0) = ln(2) = {f3_0}")
print(f"F₃(1) = ln(1 + e^(sin(1))) = {f3_1}")
print(f"The value of the integral V is ln({f3_1}) - ln({f3_0}) = {V}")
print(f"The value of 10000 * V is {final_value}")
print(f"The closest integer to 10000 * V is {closest_integer}")

# The final answer in the required format
# print(f'<<<{closest_integer}>>>')