# 1. Define the numerator and denominator for the probability of Theo NOT winning a single game.
# The probability of Theo winning is 1/6, so the probability of not winning is 5/6.
p_not_win_num = 5
p_not_win_den = 6

# 2. The problem asks for the probability of Theo's first win occurring after at least five games,
# which means Theo must not win the first four games.
num_games = 4

# 3. Calculate the numerator and denominator of the final probability.
final_numerator = p_not_win_num ** num_games
final_denominator = p_not_win_den ** num_games

# 4. Calculate the final probability as a floating-point number.
final_prob_float = final_numerator / final_denominator

# 5. Print the logic and the results.
print(f"The probability of Theo not winning a single game is {p_not_win_num}/{p_not_win_den}.")
print(f"We need to find the probability that Theo does not win in the first {num_games} consecutive games.")
print(f"This is calculated as ({p_not_win_num}/{p_not_win_den})^{num_games}.")
print(f"The final equation is: {p_not_win_num}^{num_games} / {p_not_win_den}^{num_games}")
print(f"Numerator: {p_not_win_num}^{num_games} = {final_numerator}")
print(f"Denominator: {p_not_win_den}^{num_games} = {final_denominator}")
print(f"The final probability is {final_numerator}/{final_denominator} which is approximately {final_prob_float}.")