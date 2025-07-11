import math

# Based on the sequential game theory analysis, the optimal probability of success
# for Agent C (the first mover) is p_C = 0.5.
# This value maximizes Agent C's winning probability against rational opponents B and A.

# At p_C = 0.5, the game unfolds as follows:
# 1. Agent C chooses a strategic variable w_C such that the probability of success is p_C = 0.5.
# 2. Agent B observes p_C = 0.5 and chooses w_B > w_C, with a success probability p_B approaching 1/3.
# 3. Agent A observes p_C = 0.5 and p_B approx 1/3, and chooses w_A < w_C with p_A approaching 1.
# The resulting winning probabilities for each agent are 1/3.

# C's optimal choice
p_C = 0.5

# The problem asks for floor(100 * p_C)
result = math.floor(100 * p_C)

print(f"Agent C's optimal probability choice is p_C = {p_C}")
print(f"The calculation is floor(100 * {p_C})")
print(f"The final result is: {result}")
