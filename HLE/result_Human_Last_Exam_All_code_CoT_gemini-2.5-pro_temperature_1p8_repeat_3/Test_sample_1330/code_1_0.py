import math

# --- Given parameters ---
# Number of pages for relation P
B_P = 80
# Number of pages for relation Q
B_Q = 65
# Number of available memory buffer pages
M = 15

# --- BNLJ Cost Calculation ---
# To minimize cost, the smaller relation should be the outer relation.
if B_P < B_Q:
    B_outer = B_P
    B_inner = B_Q
    outer_name = "P"
    inner_name = "Q"
else:
    B_outer = B_Q
    B_inner = B_P
    outer_name = "Q"
    inner_name = "P"
    
print(f"Choosing the smaller relation as the outer loop to minimize cost.")
print(f"Outer relation: {outer_name} ({B_outer} pages)")
print(f"Inner relation: {inner_name} ({B_inner} pages)\n")

# Number of buffer pages available for the outer relation block.
# We reserve one page for the inner relation and one for the output.
buffer_for_outer = M - 2

# Number of blocks the outer relation is partitioned into.
# Each block (except possibly the last) will be of size 'buffer_for_outer'.
num_outer_blocks = math.ceil(B_outer / buffer_for_outer)

# The total cost is the cost to scan the outer relation once, plus the cost
# of scanning the inner relation for each block of the outer relation.
total_cost = B_outer + (num_outer_blocks * B_inner)

# --- Print the results step-by-step ---
print("Block Nested Loop Join (BNLJ) Cost Calculation:")
print(f"Cost = B(outer) + (ceil(B(outer) / (M - 2)) * B(inner))")
print(f"Cost = {B_outer} + (ceil({B_outer} / ({M} - 2)) * {B_inner})")
print(f"Cost = {B_outer} + (ceil({B_outer} / {buffer_for_outer}) * {B_inner})")
print(f"Cost = {B_outer} + ({num_outer_blocks} * {B_inner})")
print(f"Cost = {B_outer} + {num_outer_blocks * B_inner}")
print(f"Total Minimum Cost = {total_cost} I/O operations")
<<<465>>>