import math

# Step 1: Define the probability of Theo winning and not winning a single game.
# As derived from the random walk analysis, the probability of Theo winning (p_T) is 1/6.
# The probability of Theo NOT winning is 1 - p_T.
prob_not_T_numerator = 5
prob_not_T_denominator = 6

# Step 2: Define the number of games Theo must not win.
# The condition "first win only after at least five games" means
# Theo does not win in games 1, 2, 3, and 4.
num_games = 4

# Step 3: Calculate the final probability.
# The probability of not winning 'num_games' in a row is (5/6)^4.
result_numerator = prob_not_T_numerator ** num_games
result_denominator = prob_not_T_denominator ** num_games

final_probability = result_numerator / result_denominator

# Step 4: Print the final equation and the result.
print(f"The final equation is ({prob_not_T_numerator}/{prob_not_T_denominator})^{num_games}")
print(f"The calculation is: {prob_not_T_numerator}^{num_games} / {prob_not_T_denominator}^{num_games} = {result_numerator} / {result_denominator}")
print(f"The final probability is: {final_probability}")
