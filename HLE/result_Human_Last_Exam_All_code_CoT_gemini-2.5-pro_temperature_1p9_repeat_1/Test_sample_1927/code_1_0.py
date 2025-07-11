# The number of nodes at each level of the Kripke frame is determined by the structure of the formula.
# Level 0: The root of the countermodel.
nodes_level_0 = 1

# Level 1: To falsify the subformula P_1, we need two distinct successor nodes.
nodes_level_1 = 2

# Level 2: Each of the two nodes at level 1 must falsify P_0, requiring two successors each.
nodes_level_2 = 2 * 2

# The total number of nodes is the sum of nodes from all levels.
total_nodes = nodes_level_0 + nodes_level_1 + nodes_level_2

# Output the equation for the total number of nodes.
print(f"The total number of nodes is the sum of nodes at each level of the Kripke frame.")
print(f"Number of nodes: {nodes_level_0} (root) + {nodes_level_1} (level 1) + {nodes_level_2} (level 2) = {total_nodes}")
