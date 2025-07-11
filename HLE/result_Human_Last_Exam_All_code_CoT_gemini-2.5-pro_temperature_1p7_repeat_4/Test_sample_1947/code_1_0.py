# Based on the step-by-step derivation, the coefficients are determined as follows:
# c1: Contribution from walks on a single edge (P_2 graph). For each edge, there are 2 walks (one starting from each endpoint). So, c1 = 2.
# c2: Contribution from walks on a triangle (K_3 graph). A tree-like walk cannot have a cycle in its set of distinct edges.
#     Walks utilizing all edges of a K_3 are not tree-like. Walks utilizing a subset of edges (1 or 2) are already
#     counted in the P_2 and P_3 cases. Thus, there is no additional contribution from k. So, c2 = 0.
# c3: Contribution from walks on a path of 3 edges (P_4 graph). Each of the 3 edges is traversed twice.
#     The derivation shows there are 6 such unique walks for each P_4 subgraph. So, c3 = 6.
# c4: Contribution from walks on a path of 2 edges (P_3 graph). The two edges are traversed (4,2) or (2,4) times.
#     The derivation shows there are 12 such unique walks for each P_3 subgraph. So, c4 = 12.
# c5: Contribution from walks on a star graph with 3 edges (K_1,3 graph). Each of the 3 edges is traversed twice.
#     The derivation shows there are 12 such unique walks for each K_1,3 subgraph. So, c5 = 12.

c1 = 2
c2 = 0
c3 = 6
c4 = 12
c5 = 12

# The problem asks for the coefficients c_1, c_2, c_3, c_4, c_5 in this order.
print(c1)
print(c2)
print(c3)
print(c4)
print(c5)