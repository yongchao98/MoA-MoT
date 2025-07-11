import math

# Calculate F3(1) and F3(0)
# F3(x) = ln(1 + e^(sin(x)))
f3_at_1 = math.log(1 + math.exp(math.sin(1)))
f3_at_0 = math.log(2)

# The integral V = ln(F3(1)) - ln(F3(0))
V = math.log(f3_at_1) - math.log(f3_at_0)

# Calculate the final result
final_result = round(10000 * V)

# Output the equation and the final result
print("The integral is V = ln(F₃(1)) - ln(F₃(0))")
print(f"F₃(1) = ln(1 + e^(sin(1))) = {f3_at_1}")
print(f"F₃(0) = ln(2) = {f3_at_0}")
print(f"V = ln({f3_at_1}) - ln({f3_at_0}) = {V}")
print(f"The value is V = {V}")
print(f"The closest integer to 10000 * V is {final_result}")
