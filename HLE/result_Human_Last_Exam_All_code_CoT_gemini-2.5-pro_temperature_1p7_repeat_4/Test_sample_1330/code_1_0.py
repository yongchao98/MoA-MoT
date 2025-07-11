import math

# Given parameters
B_P = 80  # Number of pages for relation P
B_Q = 65  # Number of pages for relation Q
M = 15    # Number of available memory buffer pages

# The block size for the outer loop is M-2
block_size = M - 2

# --- Scenario 1: P is the outer relation ---
num_outer_blocks_p = math.ceil(B_P / block_size)
cost_p_outer = B_P + (num_outer_blocks_p * B_Q)

# --- Scenario 2: Q is the outer relation ---
num_outer_blocks_q = math.ceil(B_Q / block_size)
cost_q_outer = B_Q + (num_outer_blocks_q * B_P)

# --- Determine the minimum cost and print the calculation ---
print("To calculate the minimum cost for the Block Nested Loop Join, we must decide which relation to use as the outer relation.")
print("The cost is minimized when the smaller relation is chosen as the outer relation.")
print(f"Comparing the sizes, P has {B_P} pages and Q has {B_Q} pages. Q is smaller and should be the outer relation.")
print("\nCalculating the minimum cost with Q as the outer relation:\n")
print(f"Cost = B(Q) + (ceil(B(Q) / (M - 2)) * B(P))")
print(f"Cost = {B_Q} + (ceil({B_Q} / ({M} - 2)) * {B_P})")
print(f"Cost = {B_Q} + (ceil({B_Q} / {block_size}) * {B_P})")
print(f"Cost = {B_Q} + ({num_outer_blocks_q} * {B_P})")
print(f"Cost = {B_Q} + {num_outer_blocks_q * B_P}")
print(f"Cost = {cost_q_outer}")

print("\n---")
print("For completeness, the cost with P as the outer relation would be:")
print(f"Cost = {B_P} + (ceil({B_P} / {block_size}) * {B_Q}) = {B_P} + ({num_outer_blocks_p} * {B_Q}) = {cost_p_outer}")
print(f"The minimum cost is indeed when Q is the outer relation.")

# The final answer in the specified format is derived from the minimum cost calculation.
# However, the user asked to not output the final answer tag directly in the code block.
# print(f"\n<<<{cost_q_outer}>>>")