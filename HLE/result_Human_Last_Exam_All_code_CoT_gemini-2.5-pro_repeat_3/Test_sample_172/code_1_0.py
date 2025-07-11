import math

# n is the number of vertices in the cycle
n = 100
# m is the number of parallel edges for each original edge
m = 1234

# The correspondence chromatic number is given by the formula ceil((n * m) / (n - 1))
numerator = n * m
denominator = n - 1

# Calculate the result
result = math.ceil(numerator / denominator)

# Output the equation with all the numbers
print(f"The final equation is ceil(({n} * {m}) / ({n} - 1))")
print(f"This evaluates to ceil({numerator} / {denominator})")
print(f"The correspondence chromatic number is: {result}")
