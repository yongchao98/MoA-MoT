import math

def F3(x):
  """
  Calculates the value of F3(x) = ln(1 + e^(sin(x))).
  """
  f0_x = x
  f1_x = math.sin(f0_x)
  f2_x = math.exp(f1_x)
  f3_x = math.log(1 + f2_x)
  return f3_x

# The integral of F3'(x)/F3(x) is ln(F3(x)).
# By the Fundamental Theorem of Calculus, the definite integral from 0 to 1 is:
# V = ln(F3(1)) - ln(F3(0))

# Calculate the values for the final equation
F3_at_1 = F3(1)
F3_at_0 = F3(0)

# The value of the integral V
V = math.log(F3_at_1) - math.log(F3_at_0)

# The target value is 10000 * V
target_value = 10000 * V

# The closest integer to the target value
result = round(target_value)

print("Step 1: Calculate the terms F3(1) and F3(0)")
print(f"F3(1) = ln(1 + e^sin(1)) = {F3_at_1}")
print(f"F3(0) = ln(1 + e^sin(0)) = {F3_at_0}")
print("\nStep 2: Calculate the value of the integral V = ln(F3(1)) - ln(F3(0))")
print(f"V = ln({F3_at_1}) - ln({F3_at_0}) = {V}")
print("\nStep 3: Calculate the final value, which is the closest integer to 10000 * V")
print(f"10000 * V = 10000 * {V} = {target_value}")
print(f"The closest integer is: {result}")