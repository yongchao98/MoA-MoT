import math

# Step 1: Define the function F3(x)
def f3(x):
  """Calculates the value of F3(x) = ln(1 + e^(sin(x)))"""
  return math.log(1 + math.exp(math.sin(x)))

# Step 2: Calculate the values of F3(1) and F3(0)
# These are the numbers in the final equation V = ln(F3(1)/F3(0))
f3_at_1 = f3(1)
f3_at_0 = f3(0)

print(f"The integral is ln(F3(1)/F3(0))")
print(f"Value of F3(1) = ln(1 + e^(sin(1))): {f3_at_1}")
print(f"Value of F3(0) = ln(1 + e^(sin(0))): {f3_at_0}")

# Step 3: Calculate the value of the integral V
# V = ln(F3(1) / F3(0))
V = math.log(f3_at_1 / f3_at_0)
print(f"The value of the integral V is: {V}")

# Step 4: Calculate the final result as the closest integer to 10000 * V
result = round(10000 * V)
print(f"The closest integer to 10000 * V is: {result}")
