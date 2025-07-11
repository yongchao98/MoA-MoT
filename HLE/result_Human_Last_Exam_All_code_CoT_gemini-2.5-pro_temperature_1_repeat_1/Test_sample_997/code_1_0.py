import math

# Based on the logical derivation, the optimal probability p_C that Agent C should choose is 1.
# Agent C, as the first mover, can preempt the game by choosing the maximum possible success probability.
# If C chooses p_C = 1, no other agent can choose a higher probability, so C wins.
# If C chooses any p_C < 1, the next player (Agent B) will rationally choose p_B = 1 and win.
# Therefore, C's only winning strategy is to choose p_C = 1.

p_C = 1
multiplier = 100

# The problem asks for the floor of 100 * p_C.
result = multiplier * p_C
final_answer = math.floor(result)

# As requested, printing the numbers in the final equation.
print(f"The equation is: floor({multiplier} * {p_C}) = floor({result})")
print(f"The final answer is: {final_answer}")
