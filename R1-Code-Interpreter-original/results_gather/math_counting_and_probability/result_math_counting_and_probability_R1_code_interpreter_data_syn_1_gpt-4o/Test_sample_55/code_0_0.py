import math

# Number of balls
n = 5
# Number of boxes
k = 2

# Calculate the number of ways using the stars and bars method
ways = math.comb(n + k - 1, k - 1)
print(ways)