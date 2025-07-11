import math

# Define the functions for clarity, although not strictly necessary for the final calculation
def F0(x):
    return x

def F1(x):
    return math.sin(F0(x))

def F2(x):
    return math.exp(F1(x))

def F3(x):
    return math.log(1 + F2(x))

# The integral ∫ (F₃'(x)/F₃(x)) dx from 0 to 1 evaluates to ln(F₃(1)) - ln(F₃(0))

# Calculate F₃(1) and F₃(0)
F3_at_1 = F3(1)
F3_at_0 = F3(0)

# The final equation for the integral V is: V = ln(F₃(1)) - ln(F₃(0))
# Let's print the numbers in this equation
print(f"The value of F₃(1) is: {F3_at_1}")
print(f"The value of F₃(0) is: {F3_at_0}")

# Calculate the value of the integral, V
V = math.log(F3_at_1) - math.log(F3_at_0)
print(f"The value of the integral V is: {V}")

# Calculate the final result as requested
final_result = round(10000 * V)
print(f"The closest integer to 10000 * V is: {final_result}")

print("<<<{}>>>".format(final_result))