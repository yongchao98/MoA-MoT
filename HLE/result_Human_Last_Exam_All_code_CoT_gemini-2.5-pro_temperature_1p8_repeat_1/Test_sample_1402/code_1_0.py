import math

# Step 1: Define the probability of Theo winning a single game.
# Based on the random walk model (Gambler's Ruin problem), the probability
# of Theo winning a single game (P_tw) has been calculated to be 1/6.
# This is derived from P_tw = (1/2) * P(Theo wins | H-T=1) + (1/2) * P(Theo wins | H-T=-1)
# which resolves to P_tw = (1/2) * 0 + (1/2) * (1/3) = 1/6.
p_tw_num = 1
p_tw_den = 6

# Step 2: Calculate the probability of Theo NOT winning a single game.
p_not_tw_num = p_tw_den - p_tw_num
p_not_tw_den = p_tw_den

# Step 3: Define the number of games Theo must not win consecutively.
# The question asks for the probability of Theo's first win occurring only
# after at least five games, which means he must not win the first four.
num_games = 4

# Step 4: Calculate the final probability.
# Since each game is an independent event, we raise the probability of
# Theo not winning to the power of the number of games.
final_prob_numerator = p_not_tw_num ** num_games
final_prob_denominator = p_not_tw_den ** num_games

# Step 5: Print the results and the equation components.
print("The problem asks for the probability that Theo wins for the first time only after at least five games.")
print("This is equivalent to the probability that Theo does not win any of the first four games.")
print(f"The probability of Theo winning a single game is: {p_tw_num}/{p_tw_den}")
print(f"Therefore, the probability of Theo NOT winning a single game is: {p_not_tw_num}/{p_not_tw_den}")
print(f"The number of games he must not win in a row is: {num_games}")
print("\nThe final probability is the probability of not winning, raised to the power of the number of games.")
print(f"Final Equation: ({p_not_tw_num}/{p_not_tw_den}) ^ {num_games} = {final_prob_numerator} / {final_prob_denominator}")

final_answer = final_prob_numerator / final_prob_denominator
print(f"The decimal result is: {final_answer}")

# Output the final answer in the requested format.
# The value is 625 / 1296.
final_fraction = f"{final_prob_numerator}/{final_prob_denominator}"
