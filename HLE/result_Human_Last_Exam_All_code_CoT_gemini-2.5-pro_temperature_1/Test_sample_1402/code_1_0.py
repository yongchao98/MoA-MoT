import fractions

# Step 1: Define the probability of the events in a single game.
# Based on solving the random walk model (Gambler's Ruin problem), the probabilities for a single game are:
# P(Alex wins) = 1/6
# P(Theo wins) = 1/6
# P(Draw) = 2/3
# We are interested in the probability of Theo winning.
p_theo_wins = fractions.Fraction(1, 6)

# Step 2: Calculate the probability of Theo NOT winning a single game.
# This is 1 minus the probability that he does win.
p_theo_not_wins = 1 - p_theo_wins

# Step 3: Calculate the probability for the multi-game scenario.
# The event "Theo wins for the first time only after at least five games"
# is equivalent to "Theo does not win in the first four games".
# Since each game is an independent event, we raise the probability of
# Theo not winning to the power of 4.
num_games = 4
final_probability = p_theo_not_wins ** num_games

# Step 4: Print the final equation with all the numbers.
# The equation is (P(Theo not wins))^4 = Final Probability
# which is (5/6)^4 = 625/1296.
print("The probability that Theo does not win a single game is:")
print(f"1 - {p_theo_wins.numerator}/{p_theo_wins.denominator} = {p_theo_not_wins.numerator}/{p_theo_not_wins.denominator}")
print("\nThe probability that Theo does not win in the first 4 games is calculated as:")
print(f"({p_theo_not_wins.numerator}/{p_theo_not_wins.denominator})^{num_games} = {final_probability.numerator}/{final_probability.denominator}")
print(f"\nAs a decimal, this probability is approximately {float(final_probability):.6f}")
