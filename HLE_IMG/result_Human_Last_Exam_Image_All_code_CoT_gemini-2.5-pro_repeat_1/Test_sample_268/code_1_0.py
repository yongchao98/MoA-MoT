import math
from fractions import Fraction

# The final structure is Level 3 (L=3).
# L=0 is a single node.
# L=1 is the 4-node "blue triad" (a K4 graph).
# L=2 is the 16-node "blue/green" structure.
# L=3 is the 64-node "blue/green/red" structure.
LEVEL = 3

# Step 1: Calculate the degree of the central node (k_c).
# The central node at Level L is connected to all other nodes in the graph.
# The total number of nodes at Level L is N(L) = 4^L.
# The degree of the central node, k_c, is N(L) - 1.
num_nodes = 4**LEVEL
k_c = num_nodes - 1

# Step 2: Calculate the number of edges between the neighbors of the central node (E_c).
# This is equivalent to the total number of edges in the graph minus the edges
# connected to the central node (which is k_c).
# E_c = E_total - k_c
# We need a formula for E_total(L), the total edges in a Level L graph.
# Recurrence relation: E_total(L) = 4 * E_total(L-1) + 3 * (number of nodes in L-1 graph)
# E_total(L) = 4 * E_total(L-1) + 3 * 4^(L-1)
# Base case: L=1, the graph is a K4, so E_total(1) = 6.

E_total_L_minus_1 = 6  # E_total for L=1
for l in range(2, LEVEL + 1):
    nodes_L_minus_1 = 4**(l - 1)
    E_total_L = 4 * E_total_L_minus_1 + 3 * nodes_L_minus_1
    E_total_L_minus_1 = E_total_L

E_total = E_total_L_minus_1
E_c = E_total - k_c

# Step 3: Calculate the clustering coefficient (C_c).
# Formula: C_c = (2 * E_c) / (k_c * (k_c - 1))
numerator = 2 * E_c
denominator = k_c * (k_c - 1)
C_c = numerator / denominator

print(f"Analysis of the Level-{LEVEL} ({num_nodes}-node) graph:")
print("-" * 40)

print(f"1. Degree of the central node (k_c):")
print(f"   The central node is connected to all other nodes.")
print(f"   k_c = Total Nodes - 1 = {num_nodes} - 1 = {k_c}")
print()

print(f"2. Edges between neighbors (E_c):")
print(f"   Total edges in the graph (E_total) = {E_total}")
print(f"   E_c = E_total - k_c = {E_total} - {k_c} = {E_c}")
print()

print(f"3. Clustering Coefficient (C_c):")
print(f"   Formula: C_c = (2 * E_c) / (k_c * (k_c - 1))")
print(f"   C_c = (2 * {E_c}) / ({k_c} * ({k_c} - 1))")
print(f"   C_c = {numerator} / ({k_c} * {k_c - 1})")
print(f"   C_c = {numerator} / {denominator}")
print()

final_fraction = Fraction(numerator, denominator)
print(f"Final Answer (as fraction): {final_fraction.numerator}/{final_fraction.denominator}")
print(f"Final Answer (as decimal): {C_c}")
print(f"<<<{'%.10f' % C_c}>>>")