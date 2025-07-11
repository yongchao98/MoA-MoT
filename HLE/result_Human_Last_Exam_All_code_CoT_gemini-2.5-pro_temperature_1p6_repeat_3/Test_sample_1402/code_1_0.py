import math

# Step 1: Define the probability of Theo winning a single game.
# Based on the random walk analysis, this probability is 1/6.
p_theo_wins = 1/6

# Step 2: Calculate the probability of Theo not winning a single game.
p_theo_not_wins = 1 - p_theo_wins

# The numerator and denominator for the probability of Theo not winning.
# p_theo_not_wins = 5/6
no_win_num = 5
no_win_den = 6

# Step 3: Define the number of games Theo must not win consecutively.
# The question asks for the probability of Theo winning for the first time
# after at least five games, which means he does not win the first four.
num_games = 4

# Step 4: Calculate the final probability by raising the single-game
# non-win probability to the power of the number of games.
# This is (5/6)^4.
result_num = no_win_num ** num_games
result_den = no_win_den ** num_games

# Step 5: Calculate the decimal value of the final probability.
final_probability = result_num / result_den

# Step 6: Print the explanation and the result.
print("The problem is to find the probability that Theo wins for the first time only after at least five games.")
print("This is the same as the probability that Theo does not win any of the first four games.")
print(f"The probability of Theo winning a single game is 1/6.")
print(f"Therefore, the probability of Theo NOT winning a single game is 1 - 1/6 = {no_win_num}/{no_win_den}.")
print(f"The probability of this happening for {num_games} consecutive games is ({no_win_num}/{no_win_den})^{num_games}.")
print(f"The calculation is: {no_win_num}^{num_games} / {no_win_den}^{num_games} = {result_num} / {result_den}.")
print(f"The final probability is approximately: {final_probability}")
