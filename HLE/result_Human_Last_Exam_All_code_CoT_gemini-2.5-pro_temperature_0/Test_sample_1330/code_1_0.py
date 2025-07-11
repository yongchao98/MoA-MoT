import math

# Step 1: Define the given parameters
# B(P): Number of pages for relation P
B_P = 80
# B(Q): Number of pages for relation Q
B_Q = 65
# M: Number of available memory buffer pages
M = 15

# Step 2: Determine the optimal outer and inner relations
# To minimize I/O cost in a BNLJ, the smaller relation should be the outer relation.
# This minimizes the number of times the larger, inner relation is scanned.
if B_P < B_Q:
    B_outer = B_P
    B_inner = B_Q
    outer_relation_name = "P"
    inner_relation_name = "Q"
else:
    B_outer = B_Q
    B_inner = B_P
    outer_relation_name = "Q"
    inner_relation_name = "P"

print(f"Given B(P) = {B_P} and B(Q) = {B_Q}, the smaller relation is {outer_relation_name}.")
print(f"To minimize cost, we choose {outer_relation_name} as the outer relation and {inner_relation_name} as the inner relation.")
print("-" * 20)

# Step 3: Calculate the components of the BNLJ cost formula
# The formula is: Cost = B(outer) + (ceil(B(outer) / (M - 2)) * B(inner))
# M-2 pages are available for the outer relation's block.
buffer_for_outer = M - 2

# Number of times the inner relation must be scanned (number of outer loops)
num_outer_loops = math.ceil(B_outer / buffer_for_outer)

# Total cost calculation
total_cost = B_outer + (num_outer_loops * B_inner)

# Step 4: Print the final equation with all the numbers
print("The minimum I/O cost is calculated as follows:")
print(f"Cost = B({outer_relation_name}) + (ceil(B({outer_relation_name}) / (M - 2)) * B({inner_relation_name}))")
print(f"Cost = {B_outer} + (ceil({B_outer} / ({M} - 2)) * {B_inner})")
print(f"Cost = {B_outer} + (ceil({B_outer} / {buffer_for_outer}) * {B_inner})")
print(f"Cost = {B_outer} + ({num_outer_loops} * {B_inner})")
print(f"Cost = {B_outer} + {num_outer_loops * B_inner}")
print(f"Cost = {total_cost}")

<<<465>>>