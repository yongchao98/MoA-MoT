import math

# Given parameters
B_P = 80  # Pages in relation P
B_Q = 65  # Pages in relation Q
M = 15    # Available memory buffer pages

# In Block Nested Loop Join (BNLJ), we want to choose the smaller relation
# as the outer relation to minimize the number of times we scan the inner relation.
# Let's compare the costs.

# --- Case 1: P as outer, Q as inner ---
cost_P_outer = B_P + math.ceil(B_P / (M - 1)) * B_Q

# --- Case 2: Q as outer, P as inner ---
cost_Q_outer = B_Q + math.ceil(B_Q / (M - 1)) * B_P

# The minimum cost is the smaller of the two.
# We will show the calculation for the optimal case where Q is the outer relation.

if cost_Q_outer < cost_P_outer:
    B_outer = B_Q
    B_inner = B_P
    outer_relation_name = "Q"
    inner_relation_name = "P"
else:
    B_outer = B_P
    B_inner = B_Q
    outer_relation_name = "P"
    inner_relation_name = "Q"

# Calculate the number of blocks for the outer loop buffer
buffer_for_outer = M - 1

# Calculate the number of outer loops (scans of the inner relation)
num_outer_loops = math.ceil(B_outer / buffer_for_outer)

# Calculate the total minimum cost
min_cost = B_outer + num_outer_loops * B_inner

print("To minimize I/O cost, we choose the smaller relation (Q) as the outer relation.")
print(f"The BNLJ cost formula is: B(outer) + ceil(B(outer) / (M - 1)) * B(inner)")
print(f"Calculation:")
print(f"Cost = B({outer_relation_name}) + ceil(B({outer_relation_name}) / (M - 1)) * B({inner_relation_name})")
print(f"Cost = {B_outer} + ceil({B_outer} / ({M} - 1)) * {B_inner}")
print(f"Cost = {B_outer} + ceil({B_outer} / {buffer_for_outer}) * {B_inner}")
print(f"Cost = {B_outer} + {num_outer_loops} * {B_inner}")
print(f"Cost = {B_outer} + {num_outer_loops * B_inner}")
print(f"Cost = {min_cost}")

print(f"\nThe minimum I/O cost is {min_cost}.")
<<<465>>>