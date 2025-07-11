import math

# Given values
# B(P): Number of pages for relation P
B_P = 80
# B(Q): Number of pages for relation Q
B_Q = 65
# M: Number of available memory buffer pages
M = 15

# For Block Nested Loop Join, one buffer page is used for the inner relation's block,
# one is for the output buffer, and the rest (M-2) are for the outer relation's blocks.
buffer_for_outer = M - 2

# Case 1: P is the outer relation, Q is the inner relation
# Number of chunks for the outer relation P
chunks_p_outer = math.ceil(B_P / buffer_for_outer)
# Cost = I/O for P + (Number of P chunks * I/O for Q)
cost_p_outer = B_P + (chunks_p_outer * B_Q)

# Case 2: Q is the outer relation, P is the inner relation
# Number of chunks for the outer relation Q
chunks_q_outer = math.ceil(B_Q / buffer_for_outer)
# Cost = I/O for Q + (Number of Q chunks * I/O for P)
cost_q_outer = B_Q + (chunks_q_outer * B_P)

# To minimize I/O cost, we choose the smaller of the two calculated costs.
# The selection predicates (e.g., P.x <= 15) do not affect the I/O cost of the BNLJ itself
# as they are applied after the pages are read into memory.
if cost_p_outer < cost_q_outer:
    # Print the equation for the minimum cost case (P as outer)
    print(f"Minimum cost is achieved with P as the outer relation.")
    print(f"Cost = B(P) + (ceil(B(P) / (M - 2)) * B(Q))")
    print(f"Equation: {B_P} + {chunks_p_outer} * {B_Q} = {cost_p_outer}")
else:
    # Print the equation for the minimum cost case (Q as outer)
    print(f"Minimum cost is achieved with Q as the outer relation.")
    print(f"Cost = B(Q) + (ceil(B(Q) / (M - 2)) * B(P))")
    print(f"Equation: {B_Q} + {chunks_q_outer} * {B_P} = {cost_q_outer}")
