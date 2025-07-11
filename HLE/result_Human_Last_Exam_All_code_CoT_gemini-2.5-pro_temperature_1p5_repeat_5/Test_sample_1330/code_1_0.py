import math

# Given parameters
B_P = 80  # Pages in relation P
B_Q = 65  # Pages in relation Q
M = 15    # Memory buffer pages

# The number of available buffers for the outer relation block
buffer_for_outer = M - 2

# --- Scenario 1: P as outer, Q as inner ---
num_blocks_P = math.ceil(B_P / buffer_for_outer)
cost_P_outer = B_P + num_blocks_P * B_Q

# --- Scenario 2: Q as outer, P as inner ---
num_blocks_Q = math.ceil(B_Q / buffer_for_outer)
cost_Q_outer = B_Q + num_blocks_Q * B_P

# --- Determine minimum cost and print the detailed calculation ---
print("Block Nested Loop Join (BNLJ) Cost Calculation")
print("-------------------------------------------------")
print(f"Pages in P, B(P) = {B_P}")
print(f"Pages in Q, B(Q) = {B_Q}")
print(f"Memory Pages, M = {M}")
print("")

if cost_Q_outer < cost_P_outer:
    print("Minimum cost is achieved with Q as the outer relation.")
    print("Cost = B(Q) + (ceil(B(Q) / (M - 2)) * B(P))")
    print(f"Cost = {B_Q} + (ceil({B_Q} / ({M} - 2)) * {B_P})")
    print(f"Cost = {B_Q} + (ceil({B_Q} / {buffer_for_outer}) * {B_P})")
    print(f"Cost = {B_Q} + ({num_blocks_Q} * {B_P})")
    cost_inner_scans = num_blocks_Q * B_P
    print(f"Cost = {B_Q} + {cost_inner_scans}")
    final_cost = cost_Q_outer
    print(f"Total I/O Cost = {final_cost}")
else:
    print("Minimum cost is achieved with P as the outer relation.")
    print("Cost = B(P) + (ceil(B(P) / (M - 2)) * B(Q))")
    print(f"Cost = {B_P} + (ceil({B_P} / ({M} - 2)) * {B_Q})")
    print(f"Cost = {B_P} + (ceil({B_P} / {buffer_for_outer}) * {B_Q})")
    print(f"Cost = {B_P} + ({num_blocks_P} * {B_Q})")
    cost_inner_scans = num_blocks_P * B_Q
    print(f"Cost = {B_P} + {cost_inner_scans}")
    final_cost = cost_P_outer
    print(f"Total I/O Cost = {final_cost}")

print("-------------------------------------------------")