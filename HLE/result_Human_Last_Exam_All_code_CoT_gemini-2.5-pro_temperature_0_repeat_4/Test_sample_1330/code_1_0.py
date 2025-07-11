import math

# --- Given Parameters ---
# Number of pages for relation P
B_P = 80
# Number of pages for relation Q
B_Q = 65
# Number of available memory buffer pages
M = 15

# --- BNLJ Cost Calculation ---
# To minimize I/O cost, the smaller relation is chosen as the outer relation.
# In this case, B(Q) < B(P), so Q is the outer relation.
B_outer = B_Q
B_inner = B_P

# The cost formula for BNLJ is: B(outer) + ceil(B(outer) / (M - 2)) * B(inner)
# M-2 is used because one buffer page is for reading blocks of the inner relation
# and one is for the output buffer.

# Number of chunks the outer relation is divided into
num_outer_chunks = math.ceil(B_outer / (M - 2))

# Total I/O cost
total_cost = B_outer + num_outer_chunks * B_inner

# --- Output the Calculation Steps ---
print("Calculating the minimum I/O cost for the Block Nested Loop Join.")
print(f"The smaller relation, Q (with {B_Q} pages), is chosen as the outer relation.")
print(f"The larger relation, P (with {B_P} pages), is the inner relation.")
print(f"The system has {M} memory pages available.")
print("\n--- Cost Formula ---")
print("Cost = B(outer) + ceil(B(outer) / (M - 2)) * B(inner)")
print("\n--- Calculation ---")
print(f"Cost = {B_outer} + ceil({B_outer} / ({M} - 2)) * {B_inner}")
print(f"Cost = {B_outer} + ceil({B_outer} / {M - 2}) * {B_inner}")
print(f"Cost = {B_outer} + {num_outer_chunks} * {B_inner}")
print(f"Cost = {B_outer} + {num_outer_chunks * B_inner}")
print(f"Total Minimum I/O Cost = {total_cost}")