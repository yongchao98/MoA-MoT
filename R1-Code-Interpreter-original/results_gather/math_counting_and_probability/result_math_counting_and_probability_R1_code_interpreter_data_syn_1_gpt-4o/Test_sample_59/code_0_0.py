import math

# Calculate the smallest integer n such that n^2 >= 50
n = math.ceil(math.sqrt(50))

# Calculate the largest integer m such that m^2 <= 250
m = math.floor(math.sqrt(250))

# Count the number of perfect squares between 50 and 250
count = m - n + 1

print(count)