import math

# The problem asks for the expected value of T, which we have derived to be E[T] = 3 - 2 * e^(1/3).
# In this final equation, we can identify the numerical components.
a = 3
b = 2
c = 1
d = 3

# The following print statements show each number in the final equation.
print(f"The derived equation for the expected value is E[T] = a - b * exp(c/d)")
print(f"The number for 'a' is: {a}")
print(f"The number for 'b' is: {b}")
print(f"The number for 'c' (numerator in exponent) is: {c}")
print(f"The number for 'd' (denominator in exponent) is: {d}")

# Now, we calculate the final value using the derived formula.
expected_value = a - b * math.exp(c / d)

# Print the final result.
print(f"\nThe calculated expected value of T is:")
print(expected_value)
