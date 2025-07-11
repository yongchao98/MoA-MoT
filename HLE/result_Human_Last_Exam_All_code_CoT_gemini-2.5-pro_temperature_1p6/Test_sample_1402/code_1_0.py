from fractions import Fraction

# The problem states that the probability of Theo winning a single game is 1/6.
prob_theo_wins = Fraction(1, 6)

# The probability that Theo does not win a game is 1 minus the probability that he wins.
prob_theo_not_wins = 1 - prob_theo_wins

# We are interested in the event that Theo does not win for the first 4 games.
num_games_without_win = 4

# The final probability is (probability of not winning) raised to the power of the number of games.
final_probability = prob_theo_not_wins ** num_games_without_win

# To be explicit about the numbers in the final equation, we define them here.
numerator_base = prob_theo_not_wins.numerator
denominator_base = prob_theo_not_wins.denominator
exponent = num_games_without_win
result_numerator = final_probability.numerator
result_denominator = final_probability.denominator

# Print the final equation with all its numbers.
print(f"The final equation for the probability is:")
print(f"({numerator_base}/{denominator_base})^{exponent} = {result_numerator}/{result_denominator}")