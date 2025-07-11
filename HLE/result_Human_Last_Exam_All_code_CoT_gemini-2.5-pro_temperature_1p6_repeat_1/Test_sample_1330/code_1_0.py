import math

# Given parameters
pages_p = 80
pages_q = 65
memory_pages = 15

# In Block Nested Loop Join, M-2 pages are available for the outer relation block.
# One page is for the inner relation's input buffer, and one is for the output buffer.
buffer_for_outer = memory_pages - 2

# Scenario 1: P as the outer relation, Q as the inner relation
# Calculate the number of times we need to loop through the outer relation P
outer_loops_p = math.ceil(pages_p / buffer_for_outer)
# Calculate the total cost for this scenario: B(P) + (loops * B(Q))
cost_p_outer = pages_p + (outer_loops_p * pages_q)

# Scenario 2: Q as the outer relation, P as the inner relation
# Calculate the number of times we need to loop through the outer relation Q
outer_loops_q = math.ceil(pages_q / buffer_for_outer)
# Calculate the total cost for this scenario: B(Q) + (loops * B(P))
cost_q_outer = pages_q + (outer_loops_q * pages_p)

# Determine the minimum cost and print the corresponding calculation
print("Calculating the minimum I/O cost for Block Nested Loop Join.")
print(f"Pages for P: {pages_p}")
print(f"Pages for Q: {pages_q}")
print(f"Memory Buffer Pages: {memory_pages}")
print(f"Buffer available for outer relation blocks: {memory_pages} - 2 = {buffer_for_outer}\n")

if cost_p_outer < cost_q_outer:
    print("Minimum cost is when P is the outer relation.")
    print(f"Number of outer loops = ceil({pages_p} / {buffer_for_outer}) = {outer_loops_p}")
    print(f"Cost = B(P) + (loops * B(Q))")
    print(f"Cost = {pages_p} + ({outer_loops_p} * {pages_q}) = {cost_p_outer}")
    min_cost = cost_p_outer
else:
    print("Minimum cost is when Q is the outer relation.")
    print(f"Number of outer loops = ceil({pages_q} / {buffer_for_outer}) = {outer_loops_q}")
    print(f"Cost = B(Q) + (loops * B(P))")
    print(f"Cost = {pages_q} + ({outer_loops_q} * {pages_p}) = {cost_q_outer}")
    min_cost = cost_q_outer

print(f"\nThe minimum I/O cost is {min_cost}.")
<<<465>>>