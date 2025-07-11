import math

# Based on the step-by-step analysis, the optimal strategy for Agent C is to
# choose the highest possible success probability to guarantee a win.
# The maximum value for p(w) is p(0) = 1.
p_C = 1.0

# The problem asks for the value of floor(100 * p_C).
value = 100 * p_C
result = math.floor(value)

# Output the explanation and the final equation as requested.
print("Step-by-step reasoning leads to the optimal choice for Agent C.")
print("1. Agent A wins if max(p_B, p_C) < 1.")
print("2. Knowing this, Agent B will choose p_B = 1 if p_C < 1, securing a win for B.")
print("3. Therefore, Agent C's only winning move is to choose p_C = 1 at the start.")
print(f"\nThe optimal probability for Agent C is p_C = {p_C}.")
print(f"We need to calculate floor(100 * p_C).")
print(f"The calculation is: floor({100} * {p_C}) = {result}")