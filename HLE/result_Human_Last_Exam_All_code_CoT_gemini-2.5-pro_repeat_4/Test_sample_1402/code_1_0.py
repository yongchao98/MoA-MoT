from fractions import Fraction

# Step 1: Define the probability of Theo winning a single game.
# Based on the random walk model of the game, the probability of Theo winning (P_T)
# is calculated as follows:
# P_T = P(first toss is Tail) * P(reach state -3 before 0 | starting from state -1)
# P_T = (1/2) * (1/3) = 1/6
p_theo_wins = Fraction(1, 6)

# Step 2: Formulate the probability for the multi-game scenario.
# The event "Theo wins for the first time only after at least five games"
# is equivalent to "Theo does not win in any of the first four games."
num_games_without_win = 4

# The probability of Theo NOT winning a single game is 1 - p_theo_wins.
p_theo_not_wins = 1 - p_theo_wins

# Step 3: Calculate the final probability.
# Since the games are independent, we raise the probability of Theo not winning
# to the power of the number of games.
final_prob_fraction = p_theo_not_wins ** num_games_without_win

# Print the final result in the required format.
print("This script calculates the probability that Theo wins for the first time only after at least five games.")
print("\nFirst, we determine the probability that Theo does NOT win a single game:")
print(f"P(Theo wins a game) = {p_theo_wins.numerator}/{p_theo_wins.denominator}")
print(f"P(Theo does not win a game) = 1 - {p_theo_wins.numerator}/{p_theo_wins.denominator} = {p_theo_not_wins.numerator}/{p_theo_not_wins.denominator}")

print("\nNext, we calculate the probability that Theo does not win for the first four games in a row.")
print("The final probability equation is:")
# Output each number in the final equation
print(f"({p_theo_not_wins.numerator}/{p_theo_not_wins.denominator}) ^ {num_games_without_win} = {final_prob_fraction.numerator} / {final_prob_fraction.denominator}")
print(f"\nThe decimal approximation is: {float(final_prob_fraction):.6f}")