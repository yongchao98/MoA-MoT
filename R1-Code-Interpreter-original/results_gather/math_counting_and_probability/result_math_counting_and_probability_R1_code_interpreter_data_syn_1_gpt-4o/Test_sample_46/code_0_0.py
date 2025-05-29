import math

# Number of ways to choose 1 meat from 10
ways_to_choose_meat = 10

# Number of ways to choose 2 cheeses from 9
ways_to_choose_cheese = math.comb(9, 2)

# Total number of different sandwiches
total_sandwiches = ways_to_choose_meat * ways_to_choose_cheese

print(total_sandwiches)