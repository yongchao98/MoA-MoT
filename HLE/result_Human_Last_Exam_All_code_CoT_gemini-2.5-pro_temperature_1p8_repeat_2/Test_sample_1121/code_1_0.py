# Plan:
# 1. Define the cost of a single round of the 15-to-1 magic state distillation protocol.
# 2. For a distance-3 code, logical Clifford gates are noisy. One round of distillation is sufficient.
# 3. For a distance-5 code, logical Clifford gates are much better. Two rounds are needed to match this fidelity.
# 4. The total cost is the sum of the costs for both scenarios.

# The number of physical non-Clifford gates consumed to produce one output state
# in a single round of the 15-to-1 distillation protocol.
cost_per_round = 15

# --- Scenario 1: Universal Quantum Computer on a distance-3 (d=3) surface code ---
# For a d=3 code operating near the error threshold, the logical Clifford gates are noisy.
# One round of distillation is a reasonable choice to improve T-gate fidelity without
# over-investing resources that would be nullified by the noisy Clifford operations.
num_rounds_d3 = 1
cost_d3 = cost_per_round
print(f"For a distance-3 code, we assume {num_rounds_d3} round of distillation is required.")
print(f"Number of non-Clifford gates for the d=3 implementation = {cost_d3}")
print("-" * 20)

# --- Scenario 2: Universal Quantum Computer on a distance-5 (d=5) surface code ---
# For a d=5 code, logical Clifford gates are much more reliable. To match this,
# the non-Clifford T-gate must have very high fidelity, requiring a second round of distillation.
# The input to the second round are states from the first round.
num_rounds_d5 = 2
# Cost for round 2 is 15 states from round 1. Each round 1 state costs 15 physical gates.
cost_d5 = cost_per_round * cost_per_round
print(f"For a distance-5 code, we assume {num_rounds_d5} rounds of distillation are required.")
print(f"Number of non-Clifford gates for the d=5 implementation = {cost_per_round} * {cost_per_round} = {cost_d5}")
print("-" * 20)

# --- Total Number of Gates ---
# The prompt asks for the total number required for the first scenario and then the second.
total_cost = cost_d3 + cost_d5
print("The total approximate number of non-Clifford gates required for both scenarios is the sum.")
print(f"Total gates = (Gates for d=3) + (Gates for d=5)")
print(f"Total gates = {cost_d3} + {cost_d5} = {total_cost}")
