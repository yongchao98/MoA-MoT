import math

# --- Given Parameters ---
# Pages in relation P
b_p = 80
# Pages in relation Q
b_q = 65
# Available memory buffer pages
m = 15

# --- Plan ---
# To minimize the I/O cost in a Block Nested Loop Join, we choose the smaller
# relation as the outer relation.
# Cost Formula: B(outer) + ceil(B(outer) / (M - 2)) * B(inner)
# The WHERE clause filters (P.x <= 15, Q.x > 5) are applied after pages are
# read into memory and do not change the I/O cost.

if b_p < b_q:
    b_outer = b_p
    b_inner = b_q
    outer_name = "P"
    inner_name = "Q"
else:
    b_outer = b_q
    b_inner = b_p
    outer_name = "Q"
    inner_name = "P"

print("--- Block Nested Loop Join (BNLJ) Cost Calculation ---")
print(f"To minimize cost, the smaller relation '{outer_name}' ({b_outer} pages) is chosen as the outer relation.")
print(f"The larger relation '{inner_name}' ({b_inner} pages) is the inner relation.")
print(f"System has M = {m} memory pages available.")
print("")

# --- Calculation Steps ---
# Buffer pages available for the outer relation block
buffer_for_outer = m - 2
print(f"1. Buffer pages for the outer block = M - 2 = {m} - 2 = {buffer_for_outer}")

# Number of blocks the outer relation will be split into
num_blocks_outer = math.ceil(b_outer / buffer_for_outer)
print(f"2. Number of blocks for outer relation '{outer_name}' = ceil(B({outer_name}) / {buffer_for_outer}) = ceil({b_outer} / {buffer_for_outer}) = {num_blocks_outer}")

# Total I/O cost calculation
# Cost = Read all pages of the outer relation once +
#        (Scan all pages of the inner relation for each block of the outer relation)
total_cost = b_outer + (num_blocks_outer * b_inner)
print(f"\n3. Total I/O Cost = B({outer_name}) + (Number of Blocks * B({inner_name}))")
print(f"   Final Equation: {b_outer} + ({num_blocks_outer} * {b_inner}) = {total_cost}")

print(f"\nThe minimum cost to execute the query is: {total_cost} I/O operations.")
<<<465>>>