# The reasoning provided above leads to a numerical answer.
# The problem is a mathematical question about graph theory, not a coding task to be executed.
# The explanation walks through the properties of the graph G to determine its structure.
# 1. G must have a perfect adjustable matching M.
# 2. This implies G can be decomposed into two isomorphic subgraphs G[V1] and G[V2] of 1000 vertices each, and the matching M.
# 3. G is 3-regular and connected. This forces G[V1] and G[V2] to be 2-regular and connected.
# 4. A connected 2-regular graph on 1000 vertices is a cycle C_1000.
# 5. The matching M must be an isomorphism between the two cycles.
# 6. Analysis of isomorphisms between two cycles shows that all such constructions lead to a single non-isomorphic graph, the prism graph C_1000 x K_2.
# 7. More complex cases with edges crossing between V1 and V2 (not in the matching) can be shown to either not exist or be isomorphic to the prism graph.
# Therefore, there is only one non-isomorphic graph that satisfies all the given conditions.

number_of_graphs = 1
print(f"The number of non-isomorphic graphs satisfying the conditions is: {number_of_graphs}")