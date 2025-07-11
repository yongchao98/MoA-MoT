import math

# Given parameters
B_P = 80  # Number of pages in relation P
B_Q = 65  # Number of pages in relation Q
M = 15    # Number of available memory buffer pages

# The number of buffer pages available for the outer relation block is M - 2
buffer_for_outer = M - 2

# --- Case 1: P is the outer relation ---
# Number of chunks P is divided into
num_chunks_P = math.ceil(B_P / buffer_for_outer)
# Total cost is scanning P once, plus scanning Q for each chunk of P
cost_P_outer = B_P + (num_chunks_P * B_Q)

# --- Case 2: Q is the outer relation ---
# Number of chunks Q is divided into
num_chunks_Q = math.ceil(B_Q / buffer_for_outer)
# Total cost is scanning Q once, plus scanning P for each chunk of Q
cost_Q_outer = B_Q + (num_chunks_Q * B_P)

# The minimum cost is the smaller of the two cases
min_cost = min(cost_P_outer, cost_Q_outer)

print("To find the minimum cost, we choose the smaller relation (Q) as the outer relation.")
print("The formula for the cost is: B(Q) + (ceil(B(Q) / (M - 2)) * B(P))")
print("Calculation:")
# Final equation with each number printed
print(f"{B_Q} + (ceil({B_Q} / ({M} - 2)) * {B_P}) = {min_cost}")

print(f"<<<{min_cost}>>>")