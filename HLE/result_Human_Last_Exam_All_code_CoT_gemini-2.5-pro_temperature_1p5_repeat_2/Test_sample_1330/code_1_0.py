import math

# --- Given Parameters ---
# Number of pages in relation P
B_P = 80
# Number of pages in relation Q
B_Q = 65
# Number of available memory buffer pages
M = 15

# --- BNLJ Cost Calculation ---

# The number of available buffer pages for the outer relation's block is M-2.
# One page is for the inner relation's input buffer, and one is for the output buffer.
buffer_for_outer = M - 2

# Case 1: P is the outer relation
num_blocks_p = math.ceil(B_P / buffer_for_outer)
cost_p_outer = B_P + (num_blocks_p * B_Q)

# Case 2: Q is the outer relation
num_blocks_q = math.ceil(B_Q / buffer_for_outer)
cost_q_outer = B_Q + (num_blocks_q * B_P)

# --- Determine the minimum cost and print the result ---
print("Calculating the minimum I/O cost for Block Nested Loop Join.")
print(f"Pages in P, B(P) = {B_P}")
print(f"Pages in Q, B(Q) = {B_Q}")
print(f"Memory Pages, M = {M}\n")
print(f"Buffer pages available for an outer block = M - 2 = {M} - 2 = {buffer_for_outer}\n")

if cost_p_outer <= cost_q_outer:
    min_cost = cost_p_outer
    outer_relation_name = "P"
    B_outer = B_P
    B_inner = B_Q
    num_blocks = num_blocks_p
else:
    min_cost = cost_q_outer
    outer_relation_name = "Q"
    B_outer = B_Q
    B_inner = B_P
    num_blocks = num_blocks_q

print(f"To minimize cost, we choose the smaller relation '{outer_relation_name}' as the outer relation.")
print("The cost formula is: B(outer) + (ceil(B(outer) / (M-2)) * B(inner))")
print("\nStep 1: Calculate the number of blocks for the outer relation.")
print(f"Number of blocks = ceil(B({outer_relation_name}) / (M-2)) = ceil({B_outer} / {buffer_for_outer}) = {num_blocks}")

print("\nStep 2: Calculate the total I/O cost.")
print(f"Cost = B({outer_relation_name}) + (Number of blocks * B(inner))")
final_equation = f"Cost = {B_outer} + ({num_blocks} * {B_inner}) = {min_cost}"
print(final_equation)

print(f"\nThe minimum I/O cost is: {min_cost}")
