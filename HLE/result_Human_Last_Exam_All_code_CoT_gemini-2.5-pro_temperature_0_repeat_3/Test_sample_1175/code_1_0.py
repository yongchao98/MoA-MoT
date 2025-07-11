import math

# Define the numbers in the final simplified expression
n1 = 3
d1 = 2
base = 10
ne = 10
de = 3
n2 = 37
d2 = 4

# Calculate the final result
result = (n1 / d1) * (base ** (ne / de)) + (n2 / d2)

# Print the final equation showing all its component numbers and the calculated result
print(f"({n1}/{d1}) * {base}^({ne}/{de}) + {n2}/{d2} = {result}")