import math

# Total number of ways to arrange the lamps
total_arrangements = math.comb(6, 3)

# Total number of ways to choose 3 lamps to turn on
total_turn_on = math.comb(6, 3)

# Fix the leftmost as blue and off, and rightmost as red and on
# Remaining lamps: 2 red, 2 blue
# Choose 1 red and 1 blue to turn on from the remaining 4 lamps
favorable_turn_on = math.comb(2, 1) * math.comb(2, 1)

# Calculate the probability
probability = favorable_turn_on / (total_arrangements * total_turn_on)

print(probability)