# The number of isomorphism classes of vertex-transitive graphs on n=8 vertices
# for each degree j from 0 to 7 is based on established results from graph theory census data.
# Let n_j be the number of such graphs with degree j.

# For n=8, the counts are as follows:
# n_0: 1 (The empty graph)
# n_1: 1 (A perfect matching, 4K_2)
# n_2: 2 (The 8-cycle C_8, and two 4-cycles 2C_4)
# n_3: 5 (This includes the cube graph, 2K_4, and three other connected cubic graphs)
# By the property that the complement of a j-regular vertex-transitive graph is a (7-j)-regular
# vertex-transitive graph, we have n_j = n_{7-j}.
# n_4 = n_3 = 5
# n_5 = n_2 = 2
# n_6 = n_1 = 1
# n_7 = n_0 = 1

n = [1, 1, 2, 5, 5, 2, 1, 1]

# Print the result in the specified format
print(f"[{n[0]}, {n[1]}, {n[2]}, {n[3]}, {n[4]}, {n[5]}, {n[6]}, {n[7]}]")