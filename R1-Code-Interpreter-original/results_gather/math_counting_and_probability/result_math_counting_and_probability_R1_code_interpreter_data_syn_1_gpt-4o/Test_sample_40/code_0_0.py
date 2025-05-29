import math

# Step 1: Choose 6 people out of 7
ways_to_choose = math.comb(7, 6)

# Step 2: Arrange 6 people in a circle
ways_to_arrange = math.factorial(6 - 1)

# Total number of distinct seating arrangements
total_arrangements = ways_to_choose * ways_to_arrange

print(total_arrangements)