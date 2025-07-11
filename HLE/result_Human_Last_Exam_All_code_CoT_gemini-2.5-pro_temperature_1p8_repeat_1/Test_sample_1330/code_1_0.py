import math

# Given parameters
b_p = 80  # Number of pages for relation P
b_q = 65  # Number of pages for relation Q
m = 15    # Number of available memory buffer pages

# In BNLJ, we use M-2 pages for the outer relation's block,
# one for the inner relation's page, and one for the output buffer.
buffer_for_outer = m - 2

# --- Case 1: P is the outer relation ---
# Calculate how many times we need to loop over the inner relation Q.
# This is the number of blocks we must divide P into.
num_blocks_p = math.ceil(b_p / buffer_for_outer)
# Cost = Read P once + (Scan Q for each block of P)
cost_p_outer = b_p + (num_blocks_p * b_q)

# --- Case 2: Q is the outer relation ---
# Calculate how many times we need to loop over the inner relation P.
# This is the number of blocks we must divide Q into.
num_blocks_q = math.ceil(b_q / buffer_for_outer)
# Cost = Read Q once + (Scan P for each block of Q)
cost_q_outer = b_q + (num_blocks_q * b_p)

# --- Determine the minimum cost and print the detailed calculation ---
print("To find the minimum cost, we evaluate both relations as the outer relation.\n")
print(f"Given: B(P) = {b_p}, B(Q) = {b_q}, Memory pages M = {m}")
print(f"Buffer pages for outer relation block = M - 2 = {m} - 2 = {buffer_for_outer}\n")

print("--- Calculation with P as Outer Relation ---")
print(f"Number of blocks for P = ceil(B(P) / (M - 2)) = ceil({b_p} / {buffer_for_outer}) = {num_blocks_p}")
print(f"Cost = B(P) + (Number of blocks for P * B(Q))")
# Showing the numbers in the final equation
print(f"Cost = {b_p} + ({num_blocks_p} * {b_q}) = {b_p} + {num_blocks_p * b_q} = {cost_p_outer}\n")


print("--- Calculation with Q as Outer Relation ---")
print(f"Number of blocks for Q = ceil(B(Q) / (M - 2)) = ceil({b_q} / {buffer_for_outer}) = {num_blocks_q}")
print(f"Cost = B(Q) + (Number of blocks for Q * B(P))")
# Showing the numbers in the final equation
print(f"Cost = {b_q} + ({num_blocks_q} * {b_p}) = {b_q} + {num_blocks_q * b_p} = {cost_q_outer}\n")


if cost_p_outer < cost_q_outer:
    min_cost = cost_p_outer
    print(f"Comparing the two costs ({cost_p_outer} vs {cost_q_outer}), the minimum cost is achieved when P is the outer relation.")
else:
    min_cost = cost_q_outer
    print(f"Comparing the two costs ({cost_p_outer} vs {cost_q_outer}), the minimum cost is achieved when Q is the outer relation.")

print(f"\nThe minimum I/O cost is {min_cost} pages.")
print(f"<<<{min_cost}>>>")