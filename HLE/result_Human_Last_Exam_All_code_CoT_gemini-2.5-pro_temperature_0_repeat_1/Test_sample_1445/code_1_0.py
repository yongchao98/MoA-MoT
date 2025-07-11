# 1. Define the maximum number of 0-blocks and 1-blocks for a 100-digit sequence.
# This occurs for a fully alternating sequence (e.g., 0101... or 1010...),
# which has 50 blocks of each digit.
max_n0_blocks = 50
max_n1_blocks = 50

# 2. Calculate the worst-case cost for the transformation path via the all-0s sequence.
# This is the sum of the maximum number of 1-blocks for the initial and target sequences.
cost_via_0s = max_n1_blocks + max_n1_blocks

# 3. Calculate the worst-case cost for the transformation path via the all-1s sequence.
# This is the sum of the maximum number of 0-blocks for the initial and target sequences.
cost_via_1s = max_n0_blocks + max_n0_blocks

# 4. The minimum number of operations 'n' is the minimum of these two worst-case path costs.
# This guarantees that there is always a path that is at most this long.
n = min(cost_via_0s, cost_via_1s)

# 5. Print the step-by-step calculation.
print(f"The maximum number of 0-blocks in a worst-case 100-digit sequence is {max_n0_blocks}.")
print(f"The maximum number of 1-blocks in a worst-case 100-digit sequence is {max_n1_blocks}.")
print("\nTo find the minimum number of operations for any transformation, we consider two paths:")
print("Path 1: Transform through an all-0s sequence.")
print(f"Worst-case cost = (max 1-blocks in initial) + (max 1-blocks in target) = {max_n1_blocks} + {max_n1_blocks} = {cost_via_0s}")
print("\nPath 2: Transform through an all-1s sequence.")
print(f"Worst-case cost = (max 0-blocks in initial) + (max 0-blocks in target) = {max_n0_blocks} + {max_n0_blocks} = {cost_via_1s}")
print("\nThe guaranteed number of operations is the minimum of these two worst-case costs.")
print(f"n = min({cost_via_0s}, {cost_via_1s}) = {n}")

print(f"\nFinal Answer: The minimum number of operations is {n}.")