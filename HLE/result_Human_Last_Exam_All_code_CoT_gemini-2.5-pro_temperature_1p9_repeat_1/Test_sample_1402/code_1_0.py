# The probability of Theo winning a single game is P(T) = 1/6.
# The probability of Theo not winning a single game is 1 - P(T) = 5/6.
# The event "Theo wins for the first time only after at least five games"
# is equivalent to "Theo does not win game 1 AND does not win game 2
# AND does not win game 3 AND does not win game 4".
# Since the games are independent, we can multiply their probabilities.

# Probability of Theo not winning one game
prob_theo_not_win_numerator = 5
prob_theo_not_win_denominator = 6

# Number of games Theo must not win in a row
num_games = 4

# Calculate the final probability P = (5/6)^4
final_prob_numerator = prob_theo_not_win_numerator ** num_games
final_prob_denominator = prob_theo_not_win_denominator ** num_games

# The equation is (5^4) / (6^4)
# We need to print each number in the final equation.
print(f"The probability is calculated by the equation: ({prob_theo_not_win_numerator}^{num_games}) / ({prob_theo_not_win_denominator}^{num_games})")
print(f"This evaluates to: {final_prob_numerator} / {final_prob_denominator}")

final_probability = final_prob_numerator / final_prob_denominator
print(f"The final probability is: {final_probability}")