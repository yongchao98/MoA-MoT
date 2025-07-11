import math

# Part 1: Calculate z, the smallest memory size in D.
# Based on the analysis above, the most optimized set of variables requires 172D.
# - digit result[158]: 158D
# - char i: 3D
# - char j: 3D
# - char current_size: 3D
# - char product: 3D
# - cent carry: 2D
# Total z = 158 + 3 + 3 + 3 + 3 + 2 = 172
z = 172

# Part 2: Calculate y, the first 3 digits of 100!
factorial_100 = math.factorial(100)
y = int(str(factorial_100)[:3])

# Print the final answer in the format z:y
print(f"{z}:{y}")