import math

# Step 1: Identify the given parameters.
B_P = 80  # Number of pages for relation P
B_Q = 65  # Number of pages for relation Q
M = 15    # Number of available memory buffer pages

# The number of available buffer pages for the outer relation is M-1.
# One page is reserved for reading the inner relation.
buffer_for_outer = M - 1

# Step 2: Calculate the cost for Scenario 1 (P as the outer relation).
# We calculate how many chunks of P we need to process.
# A chunk is M-1 pages of the outer relation.
chunks_P_outer = math.ceil(B_P / buffer_for_outer)
# The total I/O cost is reading P once, plus reading Q once for each chunk of P.
cost_P_outer = B_P + chunks_P_outer * B_Q

# Step 3: Calculate the cost for Scenario 2 (Q as the outer relation).
# We calculate how many chunks of Q we need to process.
chunks_Q_outer = math.ceil(B_Q / buffer_for_outer)
# The total I/O cost is reading Q once, plus reading P once for each chunk of Q.
cost_Q_outer = B_Q + chunks_Q_outer * B_P

# Step 4: Compare costs and print the details for the minimum cost option.
print("To find the minimum cost, we must evaluate both join orders:")
print("-" * 60)

# Print details for Scenario 1
print("Scenario 1: P is the outer relation, Q is the inner relation.")
print(f"Cost = B(P) + (ceil(B(P) / (M - 1)) * B(Q))")
print(f"Cost = {B_P} + (ceil({B_P} / ({M} - 1)) * {B_Q})")
print(f"Cost = {B_P} + ({int(chunks_P_outer)} * {B_Q}) = {int(cost_P_outer)} I/Os")
print("-" * 60)

# Print details for Scenario 2
print("Scenario 2: Q is the outer relation, P is the inner relation.")
print(f"Cost = B(Q) + (ceil(B(Q) / (M - 1)) * B(P))")
print(f"Cost = {B_Q} + (ceil({B_Q} / ({M} - 1)) * {B_P})")
print(f"Cost = {B_Q} + ({int(chunks_Q_outer)} * {B_P}) = {int(cost_Q_outer)} I/Os")
print("-" * 60)

# Determine the winner and print the final conclusion and equation.
if cost_Q_outer < cost_P_outer:
    print(f"\nConclusion: The minimum cost is {int(cost_Q_outer)} I/Os by choosing Q as the outer relation.")
    print("The final cost equation and calculation is:")
    print(f"Minimum Cost = {B_Q} + ceil({B_Q} / ({M} - 1)) * {B_P}")
    print(f"Minimum Cost = {B_Q} + {int(chunks_Q_outer)} * {B_P}")
    print(f"Minimum Cost = {B_Q} + {int(chunks_Q_outer * B_P)}")
    print(f"Minimum Cost = {int(cost_Q_outer)}")
else:
    print(f"\nConclusion: The minimum cost is {int(cost_P_outer)} I/Os by choosing P as the outer relation.")
    print("The final cost equation and calculation is:")
    print(f"Minimum Cost = {B_P} + ceil({B_P} / ({M} - 1)) * {B_Q}")
    print(f"Minimum Cost = {B_P} + {int(chunks_P_outer)} * {B_Q}")
    print(f"Minimum Cost = {B_P} + {int(chunks_P_outer * B_Q)}")
    print(f"Minimum Cost = {int(cost_P_outer)}")
