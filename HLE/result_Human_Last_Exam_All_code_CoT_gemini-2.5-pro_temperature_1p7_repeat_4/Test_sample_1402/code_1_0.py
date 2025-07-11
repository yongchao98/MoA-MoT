import fractions

# Step 1: Define the probabilities for a single game.
# The coin-tossing game can be modeled as a 1D random walk. Through solving
# the properties of this walk, we can find the probabilities of each outcome.
# The state is the difference d = (Number of Heads) - (Number of Tails).
# The game ends if d reaches 3 (Alex wins), -3 (Theo wins), or 0 (Draw).
# The calculated probability of Theo winning a single game, P(T), is 1/6.
p_T_numerator = 1
p_T_denominator = 6
p_T = fractions.Fraction(p_T_numerator, p_T_denominator)

# Step 2: Calculate the probability that Theo does NOT win a single game.
p_not_T = 1 - p_T

# Step 3: Calculate the final probability.
# The problem asks for the probability that Theo wins for the first time
# only after at least five games. This is equivalent to the probability that
# Theo does not win in any of the first four games.
num_games_without_win = 4
final_probability = p_not_T ** num_games_without_win

# Step 4: Print the results and the final equation.
print(f"The probability of Theo winning a single game is: {p_T}")
print(f"The probability of Theo NOT winning a single game is: 1 - {p_T} = {p_not_T}")
print(f"\nThe probability of Theo not winning for {num_games_without_win} consecutive games is calculated as ({p_not_T})^{num_games_without_win}.")
print("\nThe final equation is:")
print(f"({p_not_T.numerator}/{p_not_T.denominator}) ^ {num_games_without_win} = {final_probability.numerator} / {final_probability.denominator}")