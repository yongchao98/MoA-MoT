import math

# Total number of balls
total_balls = 11

# Number of white balls
white_balls = 5

# Total ways to draw 2 balls from 11
total_ways = math.comb(total_balls, 2)

# Ways to draw 2 white balls from 5
white_ways = math.comb(white_balls, 2)

# Probability that both balls are white
probability_white = white_ways / total_ways

print(probability_white)