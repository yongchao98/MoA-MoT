# Step 1: Define the probabilities for a single game.
# Based on solving the random walk model as described in the plan:
# The probability of Theo winning a single game (p_T) is 1/6.
# By symmetry, the probability of Alex winning is also 1/6.
# The probability of a draw is 1 - 1/6 - 1/6 = 4/6 = 2/3.

p_theo_wins = 1/6

# Step 2: Define the event of interest.
# "Theo wins for the first time only after at least five games" means
# Theo does NOT win in Game 1, AND does NOT win in Game 2,
# AND does NOT win in Game 3, AND does NOT win in Game 4.
# The outcome of Game 5 and beyond doesn't affect this condition.

# The probability of Theo NOT winning a single game is 1 - p_theo_wins.
p_theo_not_wins = 1 - p_theo_wins

# Step 3: Calculate the final probability.
# Since the games are independent, we multiply the probabilities.
# The probability is (P(Theo does not win))^4.
num_games_without_win = 4

# Numerator of the final probability fraction
numerator = 5
# Denominator of the final probability fraction
denominator = 6

final_numerator = numerator ** num_games_without_win
final_denominator = denominator ** num_games_without_win

# Calculate the decimal value
final_probability = p_theo_not_wins ** num_games_without_win

# Output the results, showing the numbers in the final equation.
print("The probability of Theo winning a single game is 1/6.")
print("The probability of Theo not winning a single game is 1 - 1/6 = 5/6.")
print(f"The probability that Theo does not win in the first {num_games_without_win} games is calculated as ({numerator}/{denominator})^{num_games_without_win}.")
print(f"This equals {final_numerator}/{final_denominator}.")
print(f"The final probability is: {final_probability}")
