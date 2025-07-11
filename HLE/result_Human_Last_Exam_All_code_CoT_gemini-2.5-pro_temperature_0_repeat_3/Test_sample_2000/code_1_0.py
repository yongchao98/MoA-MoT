# Define the 7 types of vertices based on their membership in hyperedges e1, e2, e3.
# v_1 means in e1 only, v_12 means in e1 and e2 only, etc.
v_1 = 'v_{1}'
v_2 = 'v_{2}'
v_3 = 'v_{3}'
v_12 = 'v_{12}'
v_13 = 'v_{13}'
v_23 = 'v_{23}'
v_123 = 'v_{123}'

# Define the three hyperedges of the schema hypergraph.
e1 = {v_1, v_12, v_13, v_123}
e2 = {v_2, v_12, v_23, v_123}
e3 = {v_3, v_13, v_23, v_123}

# Calculate the pairwise intersections of the hyperedges.
e1_intersect_e2 = e1.intersection(e2)
e1_intersect_e3 = e1.intersection(e3)
e2_intersect_e3 = e2.intersection(e3)

# The core adhesion set is the union of the pairwise intersections.
# For any GHTW decomposition, there must be a bag that contains this set.
core_bag = e1_intersect_e2.union(e1_intersect_e3).union(e2_intersect_e3)

# The size of this core bag gives a lower bound on the maximum bag size.
core_bag_size = len(core_bag)

# We can construct a decomposition where the maximum bag size is max(|e1|, |e2|, |e3|, |core_bag|).
# In this specific case, all these sizes are equal.
max_bag_size = max(len(e1), len(e2), len(e3), core_bag_size)

# The generalized hypertreewidth is the maximum bag size minus 1.
ghtw = max_bag_size - 1

print(f"The schema hypergraph has hyperedges of size {len(e1)}, {len(e2)}, and {len(e3)}.")
print(f"The union of pairwise intersections forms a core set of size {core_bag_size}.")
print(f"The maximum bag size in an optimal decomposition is {max_bag_size}.")
print(f"The maximum generalised hypertreewidth is {max_bag_size} - 1 = {ghtw}.")