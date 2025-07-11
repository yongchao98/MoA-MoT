from fractions import Fraction

# The problem is to find the probability that Theo wins for the first time
# only after at least five games. This means games 1, 2, 3, and 4
# are not won by Theo.

# Step 1: Define the probability of Theo winning a single game.
# Based on the analysis of the random walk described in the plan, the probability
# that Theo wins a single game (p) is 1/6.
p_theo_wins = Fraction(1, 6)

# Step 2: Define the probability of Theo NOT winning a single game.
# This is q = 1 - p.
q_theo_not_wins = 1 - p_theo_wins

# Step 3: Calculate the probability of Theo not winning the first four games.
# Since the games are independent, this probability is q^4.
num_games = 4
final_probability = q_theo_not_wins ** num_games

# Print the explanation and the result.
print("This problem asks for the probability that Theo wins for the first time on game 5 or later.")
print("This is equivalent to calculating the probability that Theo does not win any of the first four games.")
print("")

p_num = p_theo_wins.numerator
p_den = p_theo_wins.denominator
print(f"The probability of Theo winning a single game is p = {p_num}/{p_den}.")

q_num = q_theo_not_wins.numerator
q_den = q_theo_not_wins.denominator
print(f"The probability of Theo not winning a single game is q = 1 - p = {q_num}/{q_den}.")
print("")

final_num = final_probability.numerator
final_den = final_probability.denominator

print(f"The probability of Theo not winning 4 games in a row is q^4.")
print(f"Calculation: ({q_num}/{q_den})^4 = {q_num**4} / {q_den**4} = {final_num} / {final_den}")
print("")

print(f"The final probability is {final_num}/{final_den}, which is approximately {float(final_probability):.6f}.")
