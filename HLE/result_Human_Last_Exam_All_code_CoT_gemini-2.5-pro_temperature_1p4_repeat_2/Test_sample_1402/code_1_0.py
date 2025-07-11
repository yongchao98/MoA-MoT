import math

# --- Step 1: Define the probabilities for a single game ---
# Based on a random walk model, the probability of Theo winning a single game is 1/6.
p_theo_wins_num = 1
p_theo_wins_den = 6

# --- Step 2: Define the event ---
# The event is "Theo wins for the first time only after at least five games".
# This is equivalent to "Theo does not win in the first four games".
num_games = 4

# --- Step 3: Calculate the probability of the complementary event for a single game ---
# The probability of Theo NOT winning a single game.
# P(not T) = 1 - P(T) = 1 - 1/6 = 5/6
p_not_theo_win_num = p_theo_wins_den - p_theo_wins_num
p_not_theo_win_den = p_theo_wins_den

# --- Step 4: Calculate the final probability for independent games ---
# The probability is (P(not T))^4 because the games are independent.
final_prob_num = p_not_theo_win_num ** num_games
final_prob_den = p_not_theo_win_den ** num_games

# --- Step 5: Print the final equation ---
# The problem requires printing the numbers in the final equation.
print(f"The probability of Theo not winning a single game is {p_not_theo_win_num}/{p_not_theo_win_den}.")
print(f"The probability that Theo does not win in the first {num_games} games is the probability of not winning raised to the power of {num_games}.")
print("\nFinal Equation:")
print(f"({p_not_theo_win_num}/{p_not_theo_win_den})^{num_games} = {final_prob_num}/{final_prob_den}")
