import math

# --- Given Parameters ---
# Pages in relation P
bp = 80
# Pages in relation Q
bq = 65
# Pages in memory buffer
m = 15

# --- BNLJ Cost Calculation ---

# The number of available buffers for the outer relation block.
# We subtract 1 for the inner relation scan and 1 for the output buffer.
buffer_for_outer = m - 2

# Scenario 1: P is the outer relation, Q is the inner relation
# Calculate the number of outer loop iterations (number of "super-blocks" from P)
p_as_outer_loops = math.ceil(bp / buffer_for_outer)
cost_p_outer = bp + (p_as_outer_loops * bq)

# Scenario 2: Q is the outer relation, P is the inner relation
# Calculate the number of outer loop iterations (number of "super-blocks" from Q)
q_as_outer_loops = math.ceil(bq / buffer_for_outer)
cost_q_outer = bq + (q_as_outer_loops * bp)

# --- Determine the minimum cost and print the result ---
print("Calculating the minimum cost for the Block Nested Loop Join...")
print(f"Pages in P: B(P) = {bp}")
print(f"Pages in Q: B(Q) = {bq}")
print(f"Memory Pages: M = {m}")
print(f"Buffer pages for outer relation = M - 2 = {buffer_for_outer}\n")

print("Scenario 1: P is the outer relation")
print(f"Cost = B(P) + (ceil(B(P) / (M - 2)) * B(Q))")
print(f"Cost = {bp} + (ceil({bp} / {buffer_for_outer}) * {bq})")
print(f"Cost = {bp} + ({p_as_outer_loops} * {bq})")
print(f"Cost = {bp} + {p_as_outer_loops * bq} = {cost_p_outer}\n")


print("Scenario 2: Q is the outer relation")
print(f"Cost = B(Q) + (ceil(B(Q) / (M - 2)) * B(P))")
print(f"Cost = {bq} + (ceil({bq} / {buffer_for_outer}) * {bp})")
print(f"Cost = {bq} + ({q_as_outer_loops} * {bp})")
print(f"Cost = {bq} + {q_as_outer_loops * bp} = {cost_q_outer}\n")


if cost_q_outer < cost_p_outer:
    print("Minimum cost is achieved when Q is the outer relation.")
    min_cost = cost_q_outer
else:
    print("Minimum cost is achieved when P is the outer relation.")
    min_cost = cost_p_outer

print(f"\nThe minimum cost is {min_cost} I/O operations.")
