import math

# Step 1: Calculate the probabilities for a single game.
# Let's analyze the game as a random walk on the state D = H - T.
# The walk starts at D=0. The first toss must go to D=1 or D=-1.
# After the first toss, the game ends if D returns to 0 (draw),
# or if D reaches 3 (Alex win) or -3 (Theo win).

# Let p_k be the probability of Alex winning (reaching 3 before 0 or -3)
# starting from state D=k. From standard gambler's ruin analysis, we find that:
# - The probability of Alex winning starting from D=1 is p_1 = 1/3.
# - The probability of Alex winning starting from D=-1 is p_-1 = 0.
# The probability of Alex winning the game, P(A), is:
# P(A) = P(1st toss is H) * p_1 + P(1st toss is T) * p_-1
# P(A) = 0.5 * (1/3) + 0.5 * 0 = 1/6

p_alex_wins = 1/6

# By symmetry, the probability of Theo winning, P(T), is also 1/6.
p_theo_wins = 1/6

# The probability of a draw, P(D), is the remaining probability.
p_draw = 1 - p_alex_wins - p_theo_wins
# P(D) = 1 - 1/6 - 1/6 = 4/6 = 2/3

# Step 2: Calculate the probability that Theo's first win occurs after at least 5 games.
# This means that Theo does not win in game 1, game 2, game 3, AND game 4.

# The probability that Theo does NOT win a single game is:
# P(not T win) = 1 - P(T) = P(A) + P(D)
p_not_theo_wins = 1 - p_theo_wins
p_not_theo_wins_num = 5
p_not_theo_wins_den = 6

# The games are independent events. So, the probability that Theo does not win the
# first four games is P(not T win)^4.
num_games = 4
final_prob_num = p_not_theo_wins_num ** num_games
final_prob_den = p_not_theo_wins_den ** num_games

final_prob_decimal = final_prob_num / final_prob_den

# Step 3: Print the results.
print("--- Single Game Probabilities ---")
print(f"Probability Alex wins, P(A): 1/6")
print(f"Probability Theo wins, P(T): 1/6")
print(f"Probability of a Draw, P(D): 4/6 = 2/3")
print("-" * 33)
print("\n--- Main Problem Calculation ---")
print("The problem asks for the probability that Theo wins for the first time only after at least five games.")
print("This means Theo must not win in Game 1, Game 2, Game 3, and Game 4.")
print(f"The probability that Theo does not win a single game is 1 - P(T) = 1 - 1/6 = 5/6.")
print("\nSince the games are independent, we calculate (5/6)^4:")
print(f"({p_not_theo_wins_num}/{p_not_theo_wins_den})^ {num_games} = {p_not_theo_wins_num}^{num_games} / {p_not_theo_wins_den}^{num_games} = {final_prob_num} / {final_prob_den}")
print(f"\nThe final probability is {final_prob_num}/{final_prob_den}")
print(f"As a decimal, this is approximately: {final_prob_decimal:.5f}")
