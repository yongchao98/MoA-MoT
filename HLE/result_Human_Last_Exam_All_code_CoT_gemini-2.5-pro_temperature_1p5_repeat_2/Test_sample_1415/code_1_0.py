import math

# Step 1: Calculate the values of F₃(1) and F₃(0)
# F₃(x) = ln(1 + e^(sin(x)))
f3_1 = math.log(1 + math.exp(math.sin(1)))
f3_0 = math.log(2)

print(f"The equation for F₃(x) is ln(1 + e^sin(x))")
print(f"The value of F₃(1) is: {f3_1}")
print(f"The value of F₃(0) is: {f3_0}")

# Step 2: Evaluate the integral V = ln(F₃(1)) - ln(F₃(0))
V = math.log(f3_1) - math.log(f3_0)

# The numbers in the final equation
ln_f3_1 = math.log(f3_1)
ln_f3_0 = math.log(f3_0)

print(f"\nThe definite integral from 0 to 1 of F₃'(x)/F₃(x) is ln(F₃(1)) - ln(F₃(0))")
print(f"ln({f3_1}) - ln({f3_0})")
print(f"= {ln_f3_1} - ({ln_f3_0})")
print(f"= {V}")

# Step 3: Calculate 10000 * V and find the closest integer
final_value = 10000 * V
closest_integer = round(final_value)

print(f"\nThe value V is {V}")
print(f"We need to find the closest integer to 10000 * V, which is {final_value}")
print(f"The final answer is: {closest_integer}")