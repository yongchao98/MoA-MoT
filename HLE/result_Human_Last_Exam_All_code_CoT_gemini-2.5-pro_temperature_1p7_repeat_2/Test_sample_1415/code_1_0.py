import math

# Define the functions as per the problem description
def F0(x):
    return x

def F1(x):
    return math.sin(F0(x))

def F2(x):
    return math.exp(F1(x))

def F3(x):
    return math.log(1 + F2(x))

# The integral ∫ (F₃'(x)/F₃(x)) dx from 0 to 1 evaluates to ln(F₃(1)) - ln(F₃(0)).
# Let's calculate the values for F₃(1) and F₃(0).

# Calculate F₃(1)
f3_1 = F3(1)

# Calculate F₃(0)
f3_0 = F3(0) # This is math.log(2)

# Calculate the value of the integral, V
V = math.log(f3_1) - math.log(f3_0)

# Calculate the final result: the closest integer to 10000 * V
final_result = round(10000 * V)

# Print the step-by-step evaluation of the final equation
print(f"The integral V is ln(F₃(1)) - ln(F₃(0))")
print(f"F₃(1) = ln(1 + e^(sin(1))) = {f3_1}")
print(f"F₃(0) = ln(1 + e^(sin(0))) = {f3_0}")
print(f"V = ln({f3_1}) - ln({f3_0}) = {V}")
print(f"\nThe value of 10000 * V is {10000 * V}")
print(f"The closest integer to 10000 * V is: {final_result}")
