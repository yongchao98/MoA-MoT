import math

# This script calculates the expected value E[T] based on the derived formula.

# Step 1: Explain the derived formula for E[T].
print("The problem is to find the expected value of T = sum_{i=1 to m} Y_i / 3^i.")
print("Through analysis involving the law of total expectation and properties of order statistics,")
print("the expected value E[T] can be expressed as an infinite sum:")
print("E[T] = sum_{i=1 to inf} (1 / 3^i) * (i / (i+1)!)")

# Step 2: Explain the simplification of the sum.
print("\nThis sum can be simplified using the Taylor series for e^x.")
print("The sum evaluates to a closed-form expression: 3 - 2 * e^(1/3).")

# Step 3: Define the numbers in the final equation.
# The equation is result = a - b * exp(c/d)
a = 3
b = 2
c = 1
d = 3

# Step 4: Output the final equation with its numbers.
print(f"\nThe final equation is: {a} - {b} * e^({c}/{d})")

# Step 5: Calculate the numerical result.
# We use the math.exp() function for the exponential.
result = a - b * math.exp(c / d)

# Step 6: Print the final numerical value.
print("\nThe numerical value of E[T] is:")
print(result)
