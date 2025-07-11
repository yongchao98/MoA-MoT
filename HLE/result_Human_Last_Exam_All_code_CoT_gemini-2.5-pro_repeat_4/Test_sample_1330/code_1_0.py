import math

# --- Given Parameters ---
# Number of pages in relation P
B_P = 80
# Number of pages in relation Q
B_Q = 65
# Number of available memory buffer pages
M = 15

# --- BNLJ Cost Calculation ---
# To minimize the I/O cost for a Block Nested Loop Join, the smaller
# relation should be chosen as the outer relation.
# Here, B(Q) < B(P), so Q will be the outer relation and P will be the inner.

B_outer = B_Q
B_inner = B_P

# The cost formula is: B(outer) + (ceil(B(outer) / (M - 2)) * B(inner))
# M-2 pages are available for the outer relation's block.
buffer_for_outer = M - 2

# Calculate the number of blocks the outer relation is divided into.
# This determines how many times the inner relation must be scanned.
num_outer_blocks = math.ceil(B_outer / buffer_for_outer)

# Calculate the total I/O cost.
# Cost = (I/Os to read outer relation) + (I/Os to read inner relation for each outer block)
total_cost = B_outer + (num_outer_blocks * B_inner)

print("--- Minimum BNLJ I/O Cost Calculation ---")
print(f"Pages in P, B(P) = {B_P}")
print(f"Pages in Q, B(Q) = {B_Q}")
print(f"Memory Pages, M = {M}")
print("\nTo minimize cost, the smaller relation (Q) is the outer relation.")

print("\nCost Formula: B(Q) + (ceil(B(Q) / (M - 2)) * B(P))")

print("\nFinal Equation:")
# The final calculation is printed step-by-step as requested.
# The calculation shows: B_Q + num_outer_blocks * B_P = total_cost
print(f"{B_outer} + (ceil({B_outer} / ({M} - 2)) * {B_inner}) = {B_outer} + ({num_outer_blocks} * {B_inner}) = {B_outer + num_outer_blocks * B_inner} = {total_cost}")

print(f"\nThe minimum cost for the join is: {total_cost} I/O operations.")