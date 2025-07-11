import math

# Step 1: Calculate the value of F3(0)
# F3(0) = ln(1 + e^(sin(0))) = ln(1 + e^0) = ln(2)
f3_0 = math.log(2)

# Step 2: Calculate the value of F3(1)
# F3(1) = ln(1 + e^(sin(1)))
# math.sin() takes the angle in radians, which is correct here.
f3_1 = math.log(1 + math.exp(math.sin(1)))

# Step 3: Calculate the value of the definite integral, V.
# V = ln(F3(1)) - ln(F3(0))
v = math.log(f3_1) - math.log(f3_0)

# As requested, printing each number in the final equation for V
print("The definite integral V evaluates to ln(F3(1)) - ln(F3(0)).")
print(f"The value for F3(1) is: {f3_1}")
print(f"The value for F3(0) is: {f3_0}")
print(f"So, the value V = ln({f3_1}) - ln({f3_0}) is: {v}")

# Step 4: Calculate 10000 * V and round to the closest integer.
final_answer = round(10000 * v)

print(f"\nThe value of 10000 * V is {10000 * v}")
print(f"The closest integer to 10000 * V is: {final_answer}")

print(f'<<<{final_answer}>>>')