import math

# Given parameters
bp = 80  # Number of pages for relation P
bq = 65  # Number of pages for relation Q
m = 15   # Number of available memory buffer pages

# BNLJ reserves M-2 pages for the outer relation blocks
buffer_for_outer = m - 2

# Scenario 1: P is the outer relation, Q is the inner
num_chunks_p = math.ceil(bp / buffer_for_outer)
cost_p_outer = bp + num_chunks_p * bq

# Scenario 2: Q is the outer relation, P is the inner
num_chunks_q = math.ceil(bq / buffer_for_outer)
cost_q_outer = bq + num_chunks_q * bp

# The minimum cost is found by choosing the smaller relation as the outer one.
if cost_p_outer < cost_q_outer:
    min_cost = cost_p_outer
    outer_rel_pages = bp
    inner_rel_pages = bq
    num_chunks = num_chunks_p
    outer_rel_name = "P"
    inner_rel_name = "Q"
else:
    min_cost = cost_q_outer
    outer_rel_pages = bq
    inner_rel_pages = bp
    num_chunks = num_chunks_q
    outer_rel_name = "Q"
    inner_rel_name = "P"

print(f"To minimize I/O cost, we choose the smaller relation ({outer_rel_name}) as the outer relation.")
print(f"Cost formula: B({outer_rel_name}) + (ceil(B({outer_rel_name}) / (M - 2)) * B({inner_rel_name}))")
print("\nCalculating the minimum cost:")
print(f"Number of pages in outer relation {outer_rel_name}: {outer_rel_pages}")
print(f"Number of pages in inner relation {inner_rel_name}: {inner_rel_pages}")
print(f"Number of available buffer pages M: {m}")
print(f"Number of chunks for outer relation = ceil({outer_rel_pages} / ({m} - 2)) = {num_chunks}")
print(f"Final cost equation: {outer_rel_pages} + {num_chunks} * {inner_rel_pages} = {min_cost}")

print(f"\n<<< {min_cost} >>>")