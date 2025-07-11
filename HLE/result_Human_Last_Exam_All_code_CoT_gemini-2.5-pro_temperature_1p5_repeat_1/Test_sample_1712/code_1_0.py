import math

# Calculate the value of the first term: (ln(2))^2 / 2
term1 = (math.log(2)**2) / 2

# Calculate the value of the second term: 1 - pi^2 / 12
term2 = 1 - (math.pi**2) / 12

# The final result is the sum of the two terms
total_sum = term1 + term2

# Print the components of the final equation and the result
# The sum is given by the formula: (ln(2)^2)/2 + 1 - (pi^2)/12
print(f"The sum can be decomposed into two parts:")
print(f"Part 1: sum_{n=2 to inf} (-1)^n * S_{n-1}/n = (ln(2)^2)/2")
print(f"Value of Part 1 is: {term1}")
print("")
print(f"Part 2: sum_{n=2 to inf} (-1)^n / n^2 = 1 - (pi^2)/12")
print(f"Value of Part 2 is: {term2}")
print("")
print(f"The final sum is the sum of these two parts:")
print(f"Final Sum = {term1} + {term2}")
print(f"The value of the sum is: {total_sum}")