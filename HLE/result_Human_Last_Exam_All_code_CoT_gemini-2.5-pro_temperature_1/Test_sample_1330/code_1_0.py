import math

# Given parameters
B_P = 80  # Number of pages for relation P
B_Q = 65  # Number of pages for relation Q
M = 15    # Number of available memory buffer pages

# To minimize the cost of a Block Nested Loop Join, we choose the smaller relation as the outer relation.
# In this case, B(Q) = 65 is smaller than B(P) = 80.
# So, Q will be the outer relation and P will be the inner relation.

outer_relation_pages = B_Q
inner_relation_pages = B_P

# The number of available buffer pages for the outer relation block is M - 2.
# One page is for scanning the inner relation, and one is for the output buffer.
buffer_for_outer = M - 2

# Calculate the number of blocks the outer relation (Q) will be divided into.
# This is ceil(B(Q) / (M - 2)).
num_outer_blocks = math.ceil(outer_relation_pages / buffer_for_outer)

# Calculate the total I/O cost.
# Cost = (Read outer relation) + (For each outer block, read inner relation)
# Cost = B(Q) + (num_outer_blocks * B(P))
total_cost = outer_relation_pages + (num_outer_blocks * inner_relation_pages)

# Print the calculation steps
print("Calculating the minimum cost for Block Nested Loop Join:")
print(f"1. Choose the smaller relation as the outer relation to minimize I/O cost.")
print(f"   - Pages in P: {B_P}")
print(f"   - Pages in Q: {B_Q}")
print(f"   - Since {B_Q} < {B_P}, Q is the outer relation and P is the inner relation.")
print("\n2. Calculate the number of blocks for the outer relation (Q).")
print(f"   - Available buffer for outer relation = M - 2 = {M} - 2 = {buffer_for_outer} pages.")
print(f"   - Number of blocks = ceil(Pages in Q / Available buffer) = ceil({B_Q} / {buffer_for_outer}) = {num_outer_blocks}.")
print("\n3. Calculate the total I/O cost.")
print(f"   - Cost = (Pages in Q) + (Number of Q blocks * Pages in P)")
print(f"   - Final Equation: Cost = {outer_relation_pages} + ({num_outer_blocks} * {inner_relation_pages})")
print(f"   - Cost = {outer_relation_pages} + {num_outer_blocks * inner_relation_pages}")
print(f"   - Total Minimum Cost = {total_cost} I/O operations.")

print(f"\nFinal Answer Calculation:")
print(f"{outer_relation_pages} + (ceil({outer_relation_pages} / ({M} - 2)) * {inner_relation_pages}) = {total_cost}")
<<<465>>>