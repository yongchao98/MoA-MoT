import math

# Total number of points
n = 12
# Number of points to choose to form a plane
k = 3

# Calculate the number of combinations
number_of_planes = math.comb(n, k)

print(number_of_planes)