import math

# Given parameters
p_pages = 80
q_pages = 65
m_pages = 15

# Number of buffer pages available for the outer relation block
# M-2 because 1 page is for the inner relation and 1 is for the output buffer
buffer_for_outer = m_pages - 2

# --- Scenario 1: P as outer relation ---
# Number of blocks P is divided into
p_blocks = math.ceil(p_pages / buffer_for_outer)
# Cost = I/Os to read P + (Num of P blocks * I/Os to read Q)
cost_p_outer = p_pages + (p_blocks * q_pages)

# --- Scenario 2: Q as outer relation ---
# Number of blocks Q is divided into
q_blocks = math.ceil(q_pages / buffer_for_outer)
# Cost = I/Os to read Q + (Num of Q blocks * I/Os to read P)
cost_q_outer = q_pages + (q_blocks * p_pages)

# --- Determine the minimum cost and print the calculation ---
print("Block Nested Loop Join (BNLJ) Cost Calculation")
print("-------------------------------------------------")
print(f"Pages in P, B(P) = {p_pages}")
print(f"Pages in Q, B(Q) = {q_pages}")
print(f"Memory Buffer Pages, M = {m_pages}")
print("-" * 49)

if cost_p_outer < cost_q_outer:
    min_cost = cost_p_outer
    outer_relation_name = "P"
    inner_relation_name = "Q"
    outer_pages = p_pages
    inner_pages = q_pages
    outer_blocks = p_blocks
else:
    min_cost = cost_q_outer
    outer_relation_name = "Q"
    inner_relation_name = "P"
    outer_pages = q_pages
    inner_pages = p_pages
    outer_blocks = q_blocks

print(f"To minimize I/O cost, we choose the smaller relation ({outer_relation_name}) as the outer relation.")
print("\nCalculating the minimum cost:")
print(f"Cost = B({outer_relation_name}) + (ceil(B({outer_relation_name}) / (M - 2)) * B({inner_relation_name}))")
print(f"Cost = {outer_pages} + (ceil({outer_pages} / ({m_pages} - 2)) * {inner_pages})")
print(f"Cost = {outer_pages} + (ceil({outer_pages} / {buffer_for_outer}) * {inner_pages})")
print(f"Cost = {outer_pages} + ({outer_blocks} * {inner_pages})")
print(f"Cost = {outer_pages} + {outer_blocks * inner_pages}")
print(f"Minimum Cost = {min_cost} I/O operations")

print(f"\nFinal Answer: {min_cost}")
<<<465>>>