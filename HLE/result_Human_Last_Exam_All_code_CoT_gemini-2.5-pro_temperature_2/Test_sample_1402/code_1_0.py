# Step 1: Define the probability of the relevant single-game outcome.
# Based on the random walk analysis, the probability of Theo winning a single game is 1/6.
p_theo_wins = 1/6
numerator_pt = 1
denominator_pt = 6

# Step 2: Formulate the probability for the multi-game scenario.
# We want the probability that Theo wins for the first time after at least 5 games,
# which means he does not win any of the first 4 games.
num_games = 4

# The probability of Theo NOT winning a single game.
p_theo_not_wins_num = denominator_pt - numerator_pt
p_theo_not_wins_den = denominator_pt

# Step 3: Calculate the final probability.
# This is (P(Theo does not win))^4
final_prob_num = p_theo_not_wins_num ** num_games
final_prob_den = p_theo_not_wins_den ** num_games
final_prob_decimal = final_prob_num / final_prob_den

# Print the explanation and the final equation.
print("The probability of Theo winning a single game is 1/6.")
print("The probability that Theo wins for the first time after at least five games is the probability that he does not win the first four games.")
print(f"The probability of Theo not winning a single game is {p_theo_not_wins_num}/{p_theo_not_wins_den}.")
print("\nThe final equation is:")
print(f"P(not winning in {num_games} games) = ({p_theo_not_wins_num}/{p_theo_not_wins_den})^{num_games} = {final_prob_num}/{final_prob_den}")
print(f"\nThe calculated probability is {final_prob_num}/{final_prob_den}, which is approximately {final_prob_decimal:.4f}")
