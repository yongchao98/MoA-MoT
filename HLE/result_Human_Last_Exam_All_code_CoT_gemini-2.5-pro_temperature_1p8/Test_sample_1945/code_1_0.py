import math

# The condition for a coprime pair (p,q) to be good is p+q <= n+1.
# The limit of the probability Pr(n) is 3 / pi^2.

# The numbers in the final equation are 3 and 2 (from pi squared).
numerator = 3
denominator_base = math.pi
denominator_exponent = 2

# Calculate the final value
result = numerator / (denominator_base ** denominator_exponent)

# Print the components of the equation and the final result
print(f"The equation for the limit is: {numerator} / ({denominator_base}^_^{denominator_exponent})")
print(f"The numerical value of the limit is: {result}")