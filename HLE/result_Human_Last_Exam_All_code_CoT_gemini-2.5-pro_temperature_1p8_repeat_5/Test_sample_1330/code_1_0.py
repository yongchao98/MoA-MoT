import math

# Step 1: Define the given parameters.
# B(P): Number of pages for relation P
# B(Q): Number of pages for relation Q
# M: Number of available memory buffer pages
b_p = 80
b_q = 65
m = 15

# Step 2: Determine the optimal join order.
# To minimize the cost of a Block Nested Loop Join, we should use the smaller
# relation as the outer relation.
# B(Q) = 65 is smaller than B(P) = 80, so Q should be the outer relation.

# The selection predicates (e.g., P.x <= 15) do not reduce the number of page I/Os
# in this scenario because we are not using an index or materializing the
# intermediate results. The database must scan the complete original relations.

print("To calculate the minimum I/O cost for the Block Nested Loop Join, we must choose the smaller relation as the outer relation.")
print(f"Since relation Q ({b_q} pages) is smaller than relation P ({b_p} pages), we select Q as the outer relation.")
print("")

# Step 3: Calculate the minimum cost using the BNLJ formula.
# The cost formula is: B(outer) + (ceil(B(outer) / (M - 2)) * B(inner))
# Here, outer = Q and inner = P.
# We have M-2 pages available for the block from the outer relation.

# Number of available buffer pages for the outer relation block
buffer_for_outer = m - 2

# Number of times we need to read the entire inner relation (P).
# This corresponds to the number of blocks the outer relation (Q) is split into.
num_outer_scans = math.ceil(b_q / buffer_for_outer)

# Total cost calculation
total_cost = b_q + (num_outer_scans * b_p)

# Step 4: Print the detailed calculation.
# This fulfills the requirement to show each number in the final equation.
print("The cost formula with Q as the outer relation is:")
print("Cost = B(Q) + (ceil(B(Q) / (M - 2)) * B(P))")
print("")
print("Calculation steps:")
print(f"1. Reading the outer relation Q once: {b_q} I/Os")
print(f"2. Number of blocks for Q: ceil({b_q} / ({m} - 2)) = ceil({b_q} / {buffer_for_outer}) = {num_outer_scans}")
print(f"3. For each block of Q, scan the entire inner relation P: {num_outer_scans} * {b_p} = {num_outer_scans * b_p} I/Os")
print("")
print("Final equation for the minimum cost:")
print(f"Total Cost = {b_q} + {num_outer_scans} * {b_p}")
print(f"Total Cost = {b_q} + {num_outer_scans * b_p}")
print(f"Total Cost = {total_cost}")
print("")
print(f"The minimum cost to execute the query using Block Nested Loop Join is {total_cost} I/O operations.")