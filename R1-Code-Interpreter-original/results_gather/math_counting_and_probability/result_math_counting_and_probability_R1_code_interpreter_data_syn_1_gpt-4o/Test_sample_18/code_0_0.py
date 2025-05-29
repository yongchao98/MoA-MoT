import math

# Total permutations of the digits 1, 2, 3, 5, 5, 8
total_permutations = math.factorial(6) // math.factorial(2)

# Favorable permutations with last digit fixed as 5
favorable_permutations = math.factorial(5) // math.factorial(1)

# Probability calculation
probability = favorable_permutations / total_permutations

print(probability)