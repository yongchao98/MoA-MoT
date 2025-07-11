import math

# The user wants to find the value of foo(7).
# As analyzed above, foo(n) counts the number of group structures on a labeled set of n elements.
# For n=7, a prime number, this simplifies to the calculation of 7! / phi(7).

# Define the value of n.
n = 7

# Calculate the numerator, which is n!
numerator = math.factorial(n)

# Calculate the denominator, which is phi(n). For a prime n, phi(n) = n - 1.
denominator = n - 1

# Calculate the final result.
result = numerator // denominator

# Print the numbers that form the final equation, as requested.
print(f"{numerator} / {denominator} = {result}")