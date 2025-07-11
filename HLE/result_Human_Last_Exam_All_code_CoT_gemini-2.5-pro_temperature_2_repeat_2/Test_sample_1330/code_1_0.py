import math

# Given parameters
# B(P): Pages for relation P
B_P = 80
# B(Q): Pages for relation Q
B_Q = 65
# M: Available memory buffer pages
M = 15

# Calculate the number of pages available for the outer relation's block
buffer_for_outer = M - 2

# Scenario 1: P is the outer relation
# Calculate the number of outer loops (times we need to read the entire inner relation)
num_outer_loops_P = math.ceil(B_P / buffer_for_outer)
cost_P_outer = B_P + (num_outer_loops_P * B_Q)

# Scenario 2: Q is the outer relation
# Calculate the number of outer loops
num_outer_loops_Q = math.ceil(B_Q / buffer_for_outer)
cost_Q_outer = B_Q + (num_outer_loops_Q * B_P)

# Determine the minimum cost and print the corresponding calculation
print("Calculating the minimum I/O cost for Block Nested Loop Join...\n")
if cost_Q_outer < cost_P_outer:
    print("The minimum cost is achieved when Q is the outer relation.")
    print(f"Cost = B(Q) + (ceil(B(Q) / (M - 2)) * B(P))")
    print(f"Cost = {B_Q} + (ceil({B_Q} / ({M} - 2)) * {B_P})")
    print(f"Cost = {B_Q} + (ceil({B_Q / buffer_for_outer:.2f}) * {B_P})")
    print(f"Cost = {B_Q} + ({num_outer_loops_Q} * {B_P})")
    print(f"Cost = {B_Q} + {num_outer_loops_Q * B_P}")
    print(f"Minimum Cost = {cost_Q_outer}")
else:
    print("The minimum cost is achieved when P is the outer relation.")
    print(f"Cost = B(P) + (ceil(B(P) / (M - 2)) * B(Q))")
    print(f"Cost = {B_P} + (ceil({B_P} / ({M} - 2)) * {B_Q})")
    print(f"Cost = {B_P} + (ceil({B_P / buffer_for_outer:.2f}) * {B_Q})")
    print(f"Cost = {B_P} + ({num_outer_loops_P} * {B_Q})")
    print(f"Cost = {B_P} + {num_outer_loops_P * B_Q}")
    print(f"Minimum Cost = {cost_P_outer}")
