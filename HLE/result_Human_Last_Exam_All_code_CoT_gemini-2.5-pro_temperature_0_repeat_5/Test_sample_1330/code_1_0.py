import math

# --- Given Parameters ---
# Number of pages for relation P
B_P = 80
# Number of pages for relation Q
B_Q = 65
# Number of available memory buffer pages
M = 15

# --- BNLJ Cost Calculation ---

# In a Block Nested Loop Join, M-2 pages are available for the outer relation's block.
# One page is for reading the inner relation, and one is for the output buffer.
block_size = M - 2

# --- Scenario 1: P as the outer relation ---
# Cost = B(P) + ceil(B(P) / (M - 2)) * B(Q)
num_outer_blocks_P = math.ceil(B_P / block_size)
cost_P_outer = B_P + num_outer_blocks_P * B_Q

# --- Scenario 2: Q as the outer relation ---
# Cost = B(Q) + ceil(B(Q) / (M - 2)) * B(P)
num_outer_blocks_Q = math.ceil(B_Q / block_size)
cost_Q_outer = B_Q + num_outer_blocks_Q * B_P

# --- Determine the minimum cost and print the result ---
print("Calculating the minimum I/O cost for the Block Nested Loop Join.")
print(f"Given: B(P) = {B_P}, B(Q) = {B_Q}, M = {M}")
print(f"Buffer pages for outer block = M - 2 = {block_size}\n")

print(f"Cost with P as outer: {B_P} + ceil({B_P} / {block_size}) * {B_Q} = {cost_P_outer} I/Os")
print(f"Cost with Q as outer: {B_Q} + ceil({B_Q} / {block_size}) * {B_P} = {cost_Q_outer} I/Os\n")

min_cost = min(cost_P_outer, cost_Q_outer)

print("The minimum cost is achieved by choosing the smaller relation (Q) as the outer relation.")
print("Final minimum cost calculation:")
# The final equation requires printing each number
print(f"{B_Q} + ({B_Q} / {block_size}) rounded up * {B_P} = {min_cost}")
print(f"{B_Q} + {num_outer_blocks_Q} * {B_P} = {min_cost}")
print(f"{B_Q} + {num_outer_blocks_Q * B_P} = {min_cost}")
print(f"Minimum I/O Cost: {min_cost}")
<<<465>>>