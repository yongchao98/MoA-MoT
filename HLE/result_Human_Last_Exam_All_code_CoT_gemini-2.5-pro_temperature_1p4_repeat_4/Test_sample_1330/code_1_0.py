import math

# --- Given Parameters ---
# Number of pages in relation P
B_P = 80
# Number of pages in relation Q
B_Q = 65
# Number of available memory buffer pages
M = 15

# --- Plan ---
# In a Block Nested Loop Join (BNLJ), to minimize the I/O cost,
# we should choose the smaller relation as the outer relation.
# Cost formula: B(outer) + (ceil(B(outer) / (M - 2)) * B(inner))

# Determine outer and inner relations
if B_P < B_Q:
    outer_relation_name = "P"
    outer_relation_pages = B_P
    inner_relation_name = "Q"
    inner_relation_pages = B_Q
else:
    outer_relation_name = "Q"
    outer_relation_pages = B_Q
    inner_relation_name = "P"
    inner_relation_pages = B_P

# The number of pages available for the outer relation block in memory.
# We reserve one page for the inner relation scan and one for the output buffer.
buffer_for_outer = M - 2

# --- Calculation ---
# 1. Calculate the number of blocks the outer relation will be divided into.
num_blocks_outer = math.ceil(outer_relation_pages / buffer_for_outer)

# 2. Calculate the total I/O cost.
# This consists of reading the outer relation once, plus scanning the inner relation
# once for each block of the outer relation.
total_cost = outer_relation_pages + (num_blocks_outer * inner_relation_pages)

# --- Output ---
print("--- BNLJ Minimum Cost Calculation ---")
print(f"Pages in P, B(P): {B_P}")
print(f"Pages in Q, B(Q): {B_Q}")
print(f"Memory Buffer Pages, M: {M}\n")

print(f"To minimize cost, the smaller relation ({outer_relation_name}) is the outer relation.")
print(f"Outer Relation: {outer_relation_name}, Pages B({outer_relation_name}) = {outer_relation_pages}")
print(f"Inner Relation: {inner_relation_name}, Pages B({inner_relation_name}) = {inner_relation_pages}\n")

print(f"Buffer pages for an outer block = M - 2 = {M} - 2 = {buffer_for_outer}")
print(f"Number of blocks for outer relation '{outer_relation_name}' = ceil({outer_relation_pages} / {buffer_for_outer}) = {num_blocks_outer}\n")

print("--- Final Cost Equation ---")
print(f"Cost = B({outer_relation_name}) + (ceil(B({outer_relation_name}) / (M - 2)) * B({inner_relation_name}))")
print(f"Cost = {outer_relation_pages} + ({num_blocks_outer} * {inner_relation_pages})")
print(f"Cost = {outer_relation_pages} + {num_blocks_outer * inner_relation_pages}")
print(f"Minimum I/O Cost = {total_cost}")