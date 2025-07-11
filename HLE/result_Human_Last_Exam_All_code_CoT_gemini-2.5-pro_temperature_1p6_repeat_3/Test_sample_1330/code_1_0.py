import math

# Step 1: Define the given parameters
B_P = 80  # Number of pages for relation P
B_Q = 65  # Number of pages for relation Q
M = 15    # Number of available memory buffer pages

# The number of pages we can use for the outer relation's block
buffer_for_outer = M - 1

# Step 2: Calculate cost with P as the outer relation
# The number of blocks we need to read from the outer relation P
num_outer_blocks_P = math.ceil(B_P / buffer_for_outer)
# Cost = B(P) + (num_outer_blocks_P * B(Q))
cost_P_outer = B_P + (num_outer_blocks_P * B_Q)

print(f"Scenario 1: P is the outer relation")
print(f"Cost = B(P) + (ceil(B(P) / (M - 1)) * B(Q))")
print(f"Cost = {B_P} + (ceil({B_P} / ({M} - 1)) * {B_Q})")
print(f"Cost = {B_P} + ({num_outer_blocks_P} * {B_Q})")
print(f"Cost = {B_P} + {num_outer_blocks_P * B_Q}")
print(f"Total I/O cost with P as outer: {cost_P_outer}\n")


# Step 3: Calculate cost with Q as the outer relation
# The number of blocks we need to read from the outer relation Q
num_outer_blocks_Q = math.ceil(B_Q / buffer_for_outer)
# Cost = B(Q) + (num_outer_blocks_Q * B(P))
cost_Q_outer = B_Q + (num_outer_blocks_Q * B_P)

print(f"Scenario 2: Q is the outer relation")
print(f"Cost = B(Q) + (ceil(B(Q) / (M - 1)) * B(P))")
print(f"Cost = {B_Q} + (ceil({B_Q} / ({M} - 1)) * {B_P})")
print(f"Cost = {B_Q} + ({num_outer_blocks_Q} * {B_P})")
print(f"Cost = {B_Q} + {num_outer_blocks_Q * B_P}")
print(f"Total I/O cost with Q as outer: {cost_Q_outer}\n")


# Step 4: Determine the minimum cost
min_cost = min(cost_P_outer, cost_Q_outer)
print(f"Comparing the two scenarios, the minimum I/O cost is when Q is the outer relation.")
print(f"Minimum cost = {min_cost}")

# Final Answer
# <<<465>>>