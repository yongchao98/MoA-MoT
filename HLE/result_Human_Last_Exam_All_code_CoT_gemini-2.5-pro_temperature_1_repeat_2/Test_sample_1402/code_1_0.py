from fractions import Fraction

# Step 1: Define the problem structure.
# The problem asks for the probability that Theo wins for the first time after at least five games.
# This is equivalent to the probability that Theo does not win in the first four games.
# The games are independent events.
# Let p_T be the probability that Theo wins a single game.
# The desired probability is (1 - p_T)^4.
# Our main task is to calculate p_T.

# Step 2: Model a single game as a random walk.
# The state of the game is the difference d = H - T.
# The game ends if d reaches -3 (Theo wins), 3 (Alex wins), or 0 (Draw).
# Let P_T(d) be the probability of Theo winning starting from state d.
# We have the following boundary conditions:
# P_T(-3) = 1 (Theo has won)
# P_T(3) = 0 (Alex has won)
# P_T(0) = 0 (It's a draw)

# Step 3: Set up and solve the recurrence relations.
# For any non-terminal state d, the probability P_T(d) is the average of the probabilities
# of the two subsequent states: P_T(d) = 0.5 * P_T(d+1) + 0.5 * P_T(d-1).

# Let's solve for the states d = -2, -1, 1, 2.

# For d = 1 and d = 2:
# P_T(1) = 0.5 * P_T(2) + 0.5 * P_T(0) = 0.5 * P_T(2)
# P_T(2) = 0.5 * P_T(3) + 0.5 * P_T(1) = 0.5 * 0 + 0.5 * P_T(1) = 0.5 * P_T(1)
# Substituting the second equation into the first:
# P_T(1) = 0.5 * (0.5 * P_T(1)) = 0.25 * P_T(1)
# This implies (1 - 0.25) * P_T(1) = 0, so P_T(1) = 0.
p_T_1 = Fraction(0)

# For d = -1 and d = -2:
# P_T(-1) = 0.5 * P_T(0) + 0.5 * P_T(-2) = 0.5 * 0 + 0.5 * P_T(-2) = 0.5 * P_T(-2)
# P_T(-2) = 0.5 * P_T(-1) + 0.5 * P_T(-3) = 0.5 * P_T(-1) + 0.5 * 1
# Substituting the second equation into the first:
# P_T(-1) = 0.5 * (0.5 * P_T(-1) + 0.5) = 0.25 * P_T(-1) + 0.25
# (1 - 0.25) * P_T(-1) = 0.25 => 0.75 * P_T(-1) = 0.25
# P_T(-1) = 0.25 / 0.75 = 1/3
p_T_minus_1 = Fraction(1, 3)

# Step 4: Calculate the overall probability of Theo winning a game (p_T).
# The game starts at d=0. The first toss leads to d=1 (Heads) or d=-1 (Tails), each with probability 0.5.
# p_T = 0.5 * P_T(1) + 0.5 * P_T(-1)
p_T = Fraction(1, 2) * p_T_1 + Fraction(1, 2) * p_T_minus_1

# Step 5: Calculate the probability of Theo NOT winning a game.
p_not_T = 1 - p_T

# Step 6: Calculate the final probability.
# The probability that Theo does not win in the first four games is (p_not_T)^4.
num_games = 4
final_prob = p_not_T ** num_games

# Print the results in a clear, step-by-step manner.
print("Step 1: Determine the probability of Theo winning a single game (p_T).")
print(f"The probability of Theo winning, given the state d=H-T is 1, is P_T(1) = {p_T_1}.")
print(f"The probability of Theo winning, given the state d=H-T is -1, is P_T(-1) = {p_T_minus_1}.")
print(f"The game starts at d=0 and moves to d=1 or d=-1 with probability 1/2.")
print(f"Therefore, p_T = (1/2) * P_T(1) + (1/2) * P_T(-1) = (1/2) * {p_T_1} + (1/2) * {p_T_minus_1} = {p_T}.")
print("-" * 30)
print("Step 2: Calculate the probability of Theo NOT winning a single game.")
print(f"P(Theo does not win) = 1 - p_T = 1 - {p_T} = {p_not_T}.")
print("-" * 30)
print("Step 3: Calculate the probability of this happening for the first four games.")
print("The event 'Theo wins for the first time only after at least five games' means 'Theo does not win in the first 4 games'.")
print(f"The probability is (P(Theo does not win))^4 = ({p_not_T.numerator}/{p_not_T.denominator})^{num_games}.")
print(f"This evaluates to {p_not_T.numerator}^{num_games} / {p_not_T.denominator}^{num_games} = {final_prob.numerator} / {final_prob.denominator}.")
print(f"The final probability is {final_prob.numerator}/{final_prob.denominator}.")
