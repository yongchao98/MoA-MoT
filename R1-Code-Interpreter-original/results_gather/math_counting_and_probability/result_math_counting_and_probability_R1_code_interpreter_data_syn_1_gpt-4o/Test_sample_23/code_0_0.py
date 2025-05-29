import math

# Total number of digits
n = 4

# Frequencies of repeated digits
n1 = 2
n2 = 2

# Calculate the number of unique permutations
unique_permutations = math.factorial(n) // (math.factorial(n1) * math.factorial(n2))

print(unique_permutations)