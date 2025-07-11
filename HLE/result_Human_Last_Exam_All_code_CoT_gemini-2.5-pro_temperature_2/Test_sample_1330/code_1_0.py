import math

# Given values
pages_P = 80
pages_Q = 65
buffer_pages = 15

# Usable buffer for the outer relation block
block_size = buffer_pages - 2

# Case 1: P is the outer relation, Q is the inner relation
num_blocks_P = math.ceil(pages_P / block_size)
cost_P_outer = pages_P + (num_blocks_P * pages_Q)

# Case 2: Q is the outer relation, P is the inner relation
num_blocks_Q = math.ceil(pages_Q / block_size)
cost_Q_outer = pages_Q + (num_blocks_Q * pages_P)

# Determine the minimum cost and print the detailed calculation for that case.
print("To find the minimum cost, we evaluate both join orders.\n")
print(f"The system has {buffer_pages} buffer pages. We use {block_size} pages for the outer relation block.")

if cost_Q_outer < cost_P_outer:
    print("The minimum cost is achieved by using Q as the outer relation.\n")
    print("Calculation:")
    print(f"Cost = Pages_Q + ceil(Pages_Q / (Buffer_Pages - 2)) * Pages_P")
    print(f"Cost = {pages_Q} + ceil({pages_Q} / ({buffer_pages} - 2)) * {pages_P}")
    print(f"Cost = {pages_Q} + ceil({pages_Q} / {block_size}) * {pages_P}")
    print(f"Cost = {pages_Q} + {num_blocks_Q} * {pages_P}")
    print(f"Cost = {pages_Q} + {num_blocks_Q * pages_P}")
    print(f"Cost = {cost_Q_outer}")
    min_cost = cost_Q_outer
else:
    print("The minimum cost is achieved by using P as the outer relation.\n")
    print("Calculation:")
    print(f"Cost = Pages_P + ceil(Pages_P / (Buffer_Pages - 2)) * Pages_Q")
    print(f"Cost = {pages_P} + ceil({pages_P} / ({buffer_pages} - 2)) * {pages_Q}")
    print(f"Cost = {pages_P} + ceil({pages_P} / {block_size}) * {pages_Q}")
    print(f"Cost = {pages_P} + {num_blocks_P} * {pages_Q}")
    print(f"Cost = {pages_P} + {num_blocks_P * pages_Q}")
    print(f"Cost = {cost_P_outer}")
    min_cost = cost_P_outer

print(f"\nThe minimum I/O cost to execute the query is {min_cost}.")