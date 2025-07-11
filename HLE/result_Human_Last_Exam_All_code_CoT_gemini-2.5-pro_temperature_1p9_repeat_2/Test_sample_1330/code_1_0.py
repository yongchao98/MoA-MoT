import math

# Given parameters
b_p = 80  # Number of pages for relation P
b_q = 65  # Number of pages for relation Q
m = 15    # Number of available memory buffer pages

# The number of buffer pages available for the outer relation's block
# We reserve one page for the inner relation and one page for the output buffer.
buffer_for_outer = m - 2

# --- Scenario 1: P is the outer relation, Q is the inner relation ---
# Number of chunks for the outer relation P
chunks_p = math.ceil(b_p / buffer_for_outer)
cost_p_outer = b_p + (chunks_p * b_q)

# --- Scenario 2: Q is the outer relation, P is the inner relation ---
# Number of chunks for the outer relation Q
chunks_q = math.ceil(b_q / buffer_for_outer)
cost_q_outer = b_q + (chunks_q * b_p)

# The minimum cost is the optimal one.
min_cost = min(cost_p_outer, cost_q_outer)

# Output the explanation and the final answer
print("Block Nested Loop Join (BNLJ) Cost Calculation:")
print(f"Pages for P, B(P) = {b_p}")
print(f"Pages for Q, B(Q) = {b_q}")
print(f"Memory Pages, M = {m}")
print(f"Buffer pages available for outer relation block = M - 2 = {buffer_for_outer}")
print("\n---")
print("Scenario 1: P as outer relation")
print(f"Cost = B(P) + (ceil(B(P) / (M-2)) * B(Q))")
print(f"Cost = {b_p} + (ceil({b_p} / {buffer_for_outer}) * {b_q})")
print(f"Cost = {b_p} + ({chunks_p} * {b_q}) = {cost_p_outer}")
print("\n---")
print("Scenario 2: Q as outer relation (Optimal)")
print(f"Cost = B(Q) + (ceil(B(Q) / (M-2)) * B(P))")
print(f"Cost = {b_q} + (ceil({b_q} / {buffer_for_outer}) * {b_p})")
print(f"Cost = {b_q} + ({chunks_q} * {b_p}) = {cost_q_outer}")
print("\n---")
print(f"The minimum cost is achieved when the smaller relation (Q) is the outer relation.")
print(f"Final Minimum Cost = {min_cost} I/Os")

# Final answer in the specified format
# The calculation for the final answer is: 65 + (ceil(65 / 13) * 80) = 65 + (5 * 80) = 65 + 400 = 465
print(f"\nThe final calculation is: {b_q} + ({chunks_q} * {b_p}) = {min_cost}")
<<<465>>>