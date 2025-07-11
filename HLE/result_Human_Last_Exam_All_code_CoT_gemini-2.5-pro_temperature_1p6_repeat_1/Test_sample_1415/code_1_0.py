import math

def F0(x):
  return x

def F1(x):
  return math.sin(F0(x))

def F2(x):
  return math.exp(F1(x))

def F3(x):
  return math.log(1 + F2(x))

# The integral of F3'(x)/F3(x) is ln(F3(x)).
# By the Fundamental Theorem of Calculus, the definite integral from 0 to 1 is:
# V = ln(F3(1)) - ln(F3(0))

# Calculate F3(1) and F3(0)
f3_at_1 = F3(1)
f3_at_0 = F3(0)

# The integral equation with the calculated numbers
print(f"The integral V is ln(F3(1)) - ln(F3(0))")
print(f"F3(1) = ln(1 + e^(sin(1))) = {f3_at_1}")
print(f"F3(0) = ln(1 + e^(sin(0))) = {f3_at_0}")
print(f"So, the equation is: ln({f3_at_1}) - ln({f3_at_0})")


# Calculate the value of V
V = math.log(f3_at_1) - math.log(f3_at_0)
# Note: This can also be calculated as math.log(f3_at_1 / f3_at_0)

# Calculate the final result
final_value = 10000 * V
closest_integer = round(final_value)

print(f"\nThe value of the integral V is: {V}")
print(f"The value of 10000 * V is: {final_value}")
print(f"The closest integer to 10000*V is: {closest_integer}")