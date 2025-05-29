# Number of red and white marbles
red_marbles = 3
white_marbles = 5

# Total number of marbles
total_marbles = red_marbles + white_marbles

# Probability of drawing a red marble first
prob_red_first = red_marbles / total_marbles

# Probability of drawing a white marble second
# After drawing a red marble, one less marble in total and one less red marble
prob_white_second = white_marbles / (total_marbles - 1)

# Combined probability
combined_probability = prob_red_first * prob_white_second

print(combined_probability)