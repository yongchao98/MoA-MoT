import math

# The problem is a sequential game that can be modeled as an all-pay auction.
# The analysis is done via backward induction.
#
# Let p_C, p_B, p_A be the probabilities chosen by agents C, B, and A.
# The payoff functions are:
# P(win_H) = p_H * (1 - p_M) * (1 - p_L)
# P(win_M) = p_M * (1 - p_L)
# P(win_L) = p_L
# where p_H > p_M > p_L are the chosen probabilities in decreasing order.

# Our step-by-step analysis shows:
# 1. Agent A (last mover) observes p_B and p_C and chooses p_A to maximize their win probability.
# 2. Agent B observes p_C and chooses p_B to maximize their win probability, anticipating A's move.
# 3. Agent C chooses p_C first to maximize their win probability, anticipating B's and A's moves.

# The core of the derivation is as follows:
# - If C chooses p_C > 0.5, B's optimal response is to choose p_B = 0.5.
# - This is because B's payoff of 0.25 (from choosing p_B=0.5) is greater than the payoff B would get by choosing p_B > p_C, which would be (1-p_C)^2 < 0.25.
# - Given these choices (p_C > 0.5, p_B = 0.5), A's best move is to choose p_A slightly less than 0.5.
# - The resulting order of probabilities is p_C > p_B > p_A.
# - In this case, Agent C's win probability is P(win_C) = p_C * (1 - p_B) * (1 - p_A).
# - With p_B = 0.5 and p_A approaching 0.5, this becomes P(win_C) = p_C * (1 - 0.5) * (1 - 0.5) = p_C / 4.

# To maximize this payoff, p_C/4, Agent C must choose the largest possible p_C.
# However, p_C cannot be exactly 1. If p_C = 1, C's payoff becomes 0.
# So, C must choose p_C arbitrarily close to 1 (i.e., p_C = 1 - epsilon, for a small epsilon > 0).

p_C_approaches = 1.0

# The problem asks for floor(100 * p_C).
# Since p_C must be strictly less than 1, 100 * p_C is strictly less than 100.
# For any p_C = 1 - epsilon where epsilon is small enough, 100*p_C will be between 99 and 100.
# The floor of this value is 99.
p_C_optimal_value_for_floor = 0.99999 # A value arbitrarily close to 1 but less than 1.

result = math.floor(100 * p_C_optimal_value_for_floor)

# The final equation is floor(100 * p_C)
# where p_C is the optimal probability C should choose.
# As p_C approaches 1, the value of the expression becomes 99.
hundred = 100
final_p_C = "approaches 1"

# We demonstrate the calculation with a value of p_C very close to 1.
p_C_example = 0.999999
calculation = math.floor(hundred * p_C_example)

print(f"Agent C's optimal strategy is to choose a probability p_C that is arbitrarily close to 1.")
print(f"Let's represent this as p_C -> 1.")
print(f"We need to find floor(100 * p_C).")
print(f"As p_C approaches 1 from below, 100 * p_C approaches 100 from below.")
print(f"For any value of p_C in the interval (0.99, 1), the floor of 100 * p_C will be 99.")
print(f"For example, if p_C = {p_C_example}:")
print(f"floor({hundred} * {p_C_example}) = floor({hundred * p_C_example}) = {calculation}")
print(f"Thus, the result is {result}.")
