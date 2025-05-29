import math

# Total ways to choose 2 members from 10
total_ways = math.comb(10, 2)

# Ways to choose 2 girls from 5
ways_to_choose_girls = math.comb(5, 2)

# Probability that both chosen are girls
probability = ways_to_choose_girls / total_ways

print(probability)