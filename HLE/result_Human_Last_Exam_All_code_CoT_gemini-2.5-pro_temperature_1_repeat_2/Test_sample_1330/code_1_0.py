import math

# Given parameters
B_P = 80  # Number of pages for relation P
B_Q = 65  # Number of pages for relation Q
M = 15    # Number of available memory buffer pages

# Effective number of buffer pages for the outer relation block
# We reserve one page for the inner relation and one for the output buffer.
M_effective = M - 2

# --- Scenario 1: P is the outer relation ---
# Number of passes required over the inner relation Q
passes_p_outer = math.ceil(B_P / M_effective)
# Cost = Read P once + (Number of passes * Read Q for each pass)
cost_p_outer = B_P + passes_p_outer * B_Q

# --- Scenario 2: Q is the outer relation ---
# Number of passes required over the inner relation P
passes_q_outer = math.ceil(B_Q / M_effective)
# Cost = Read Q once + (Number of passes * Read P for each pass)
cost_q_outer = B_Q + passes_q_outer * B_P

# --- Determine the minimum cost ---
if cost_q_outer < cost_p_outer:
    min_cost = cost_q_outer
    outer_relation_name = "Q"
    inner_relation_name = "P"
    B_outer = B_Q
    B_inner = B_P
    passes = passes_q_outer
else:
    min_cost = cost_p_outer
    outer_relation_name = "P"
    inner_relation_name = "Q"
    B_outer = B_P
    B_inner = B_Q
    passes = passes_p_outer

print("Block Nested Loop Join (BNLJ) Cost Calculation")
print("-------------------------------------------------")
print(f"Pages in P, B(P): {B_P}")
print(f"Pages in Q, B(Q): {B_Q}")
print(f"Memory Buffer Pages, M: {M}")
print("\nTo minimize I/O cost, we choose the smaller relation as the outer relation.")
print(f"Since B(Q) = {B_Q} is less than B(P) = {B_P}, we choose Q as the outer relation.")
print("\nCost Formula: B(outer) + ceil(B(outer) / (M - 2)) * B(inner)")
print("\nCalculating the minimum cost with Q as the outer relation:")
print(f"Cost = B(Q) + ceil(B(Q) / (M - 2)) * B(P)")
print(f"Cost = {B_outer} + ceil({B_outer} / ({M} - 2)) * {B_inner}")
print(f"Cost = {B_outer} + ceil({B_outer} / {M_effective}) * {B_inner}")
print(f"Cost = {B_outer} + {passes} * {B_inner}")
print(f"Cost = {B_outer} + {passes * B_inner}")
print(f"Total Minimum Cost = {min_cost} I/O operations")

print("\n<<<465>>>")