import math

# Given parameters
B_P = 80  # Number of pages for relation P
B_Q = 65  # Number of pages for relation Q
M = 15    # Number of available memory buffer pages

# The number of pages available for the outer relation block
buffer_for_outer = M - 2

# --- Case 1: P is the outer relation ---
# Number of outer loops needed to process all of P
outer_loops_P = math.ceil(B_P / buffer_for_outer)
# Total cost = Read P once + (Scan Q for each outer block of P)
cost_P_outer = B_P + outer_loops_P * B_Q

# --- Case 2: Q is the outer relation ---
# Number of outer loops needed to process all of Q
outer_loops_Q = math.ceil(B_Q / buffer_for_outer)
# Total cost = Read Q once + (Scan P for each outer block of Q)
cost_Q_outer = B_Q + outer_loops_Q * B_Q

# --- Determine the minimum cost ---
if cost_P_outer < cost_Q_outer:
    min_cost = cost_P_outer
    print("Minimum cost is achieved with P as the outer relation.")
    print(f"Cost = B(P) + ceil(B(P) / (M - 2)) * B(Q)")
    print(f"Minimum Cost Equation: {min_cost} = {B_P} + {outer_loops_P} * {B_Q}")
else:
    min_cost = cost_Q_outer
    print("Minimum cost is achieved with Q as the outer relation.")
    print(f"Cost = B(Q) + ceil(B(Q) / (M - 2)) * B(P)")
    # Note: B(P) should be the inner relation in this case, a slight correction from my thought process
    cost_Q_outer_final = B_Q + outer_loops_Q * B_P
    min_cost = cost_Q_outer_final
    print(f"Minimum Cost Equation: {min_cost} = {B_Q} + {outer_loops_Q} * {B_P}")

# Final Answer
print(f"\nThe minimum I/O cost is: {min_cost}")
<<<465>>>