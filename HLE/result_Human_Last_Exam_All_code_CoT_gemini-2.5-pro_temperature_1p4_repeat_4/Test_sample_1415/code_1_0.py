import math

# Define the functions to calculate F3(x)
def F0(x):
    return x

def F1(x):
    return math.sin(F0(x))

def F2(x):
    return math.exp(F1(x))

def F3(x):
    return math.log(1 + F2(x))

# The limits of the integration
x_upper = 1
x_lower = 0

# Calculate F3 at the upper and lower limits
f3_at_1 = F3(x_upper)
f3_at_0 = F3(x_lower)

# The value of the definite integral V is ln(F3(1)) - ln(F3(0))
# This can also be written as ln(F3(1) / F3(0))
V = math.log(f3_at_1) - math.log(f3_at_0)

# Calculate the final value as requested
final_value = 10000 * V
closest_integer = round(final_value)

# Print the final equation with the computed numbers
print(f"The integral V = ln(F₃(1)) - ln(F₃(0))")
print(f"V = ln({f3_at_1}) - ln({f3_at_0})")
print(f"V = {V}")
print(f"The value to report is the closest integer to 10000 * V.")
print(f"10000 * {V} = {final_value}")
print(f"The closest integer is: {closest_integer}")

print("<<<{}>>>".format(closest_integer))