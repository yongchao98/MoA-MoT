import math

# The problem is to find the constant C in the relation k_{k,inf} <= C * V.
# We found this constant to be C = 4 / pi.
# The following code calculates this value and prints the numbers in the equation.

numerator = 4
denominator_symbol = "pi"
denominator_value = math.pi

upper_bound = numerator / denominator_value

print("The upper bound for the ratio (let's call it C) between the maximum norm (k_{k,inf}) and the covolume (V) is given by the equation:")
print(f"C = {numerator} / {denominator_symbol}")
print("\nThe numbers in this final equation are:")
print(f"The numerator is: {numerator}")
print(f"The denominator, pi, is approximately: {denominator_value}")
print(f"The value of the upper bound C is approximately: {upper_bound}")
