import math

# --- Given Parameters ---
# Number of pages for relation P
B_P = 80
# Number of pages for relation Q
B_Q = 65
# Number of available memory buffer pages
M = 15

# --- BNLJ Cost Calculation ---
# The cost of a Block Nested Loop Join depends on which relation is chosen as the outer relation.
# The general formula is: Cost = B(outer) + (ceil(B(outer) / (M - 2)) * B(inner))
# We calculate the cost for both scenarios and choose the minimum.
# Note: The selection predicates (P.x <= 15, etc.) don't reduce the I/O cost here because
# the problem states "without materialization", so we must scan all original pages.

# --- Scenario 1: P is the outer relation ---
cost_P_outer = B_P + (math.ceil(B_P / (M - 2)) * B_Q)

# --- Scenario 2: Q is the outer relation ---
cost_Q_outer = B_Q + (math.ceil(B_Q / (M - 2)) * B_P)

# --- Determine the minimum cost and print the optimal strategy ---
print("Calculating the minimum I/O cost for Block Nested Loop Join (BNLJ)...")
print(f"Pages in P: {B_P}")
print(f"Pages in Q: {B_Q}")
print(f"Memory Buffer Pages: {M}\n")

print("Option 1: P as outer relation, Q as inner relation.")
print(f"Cost = B(P) + (ceil(B(P) / (M - 2)) * B(Q))")
print(f"Cost = {B_P} + (ceil({B_P} / ({M} - 2)) * {B_Q})")
print(f"Cost = {B_P} + ({math.ceil(B_P / (M - 2))} * {B_Q}) = {int(cost_P_outer)}\n")

print("Option 2: Q as outer relation, P as inner relation.")
print(f"Cost = B(Q) + (ceil(B(Q) / (M - 2)) * B(P))")
print(f"Cost = {B_Q} + (ceil({B_Q} / ({M} - 2)) * {B_P})")
passes_Q_outer = math.ceil(B_Q / (M - 2))
cost_expr_Q_outer = f"{B_Q} + ({passes_Q_outer} * {B_P})"
print(f"Cost = {cost_expr_Q_outer} = {int(cost_Q_outer)}\n")

min_cost = min(cost_P_outer, cost_Q_outer)

print("--- Conclusion ---")
print("The optimal strategy is to use the smaller relation (Q) as the outer relation.")
print(f"The minimum cost is the result of the second calculation:")
print(f"Final Equation: {B_Q} + {passes_Q_outer} * {B_P} = {int(min_cost)}")
print(f"Minimum I/O Cost: {int(min_cost)}")

<<<465>>>