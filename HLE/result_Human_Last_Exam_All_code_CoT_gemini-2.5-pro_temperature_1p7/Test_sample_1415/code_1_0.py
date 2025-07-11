import math

def F0(x):
    """Calculates F0(x) = x"""
    return x

def F1(x):
    """Calculates F1(x) = sin(F0(x))"""
    return math.sin(F0(x))

def F2(x):
    """Calculates F2(x) = e^(F1(x))"""
    return math.exp(F1(x))

def F3(x):
    """Calculates F3(x) = ln(1 + F2(x))"""
    return math.log(1 + F2(x))

# The integral of F3'(x)/F3(x) dx from 0 to 1 is V = ln(F3(1)) - ln(F3(0)).
# We need to find the values for the final equation.

# Calculate F3(0) and F3(1)
f3_at_0 = F3(0)
f3_at_1 = F3(1)

# The final equation is V = ln(f3_at_1) - ln(f3_at_0)
# We print each number in this final equation.
print("The definite integral V is calculated using the equation: ln(F3(1)) - ln(F3(0))")
print(f"Value of F3(1) = ln(1 + e^(sin(1))): {f3_at_1}")
print(f"Value of F3(0) = ln(2): {f3_at_0}")

# Calculate the value of the integral, V
V = math.log(f3_at_1) - math.log(f3_at_0)
print(f"The value of the integral V is: {V}")

# Calculate 10000 * V
final_value = 10000 * V
print(f"The value of 10000 * V is: {final_value}")

# Find the closest integer
closest_integer = round(final_value)
print(f"The final answer (closest integer to 10000 * V) is: {closest_integer}")

print("<<<{}>>>".format(closest_integer))