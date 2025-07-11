# This code is intended to be run in a SageMath environment.
from sage.all import graphs

# The number of vertices for the graphs.
N = 8

# Initialize a list to store the counts n_j for degrees j = 0 to 7.
n_j = [0] * N

# SageMath's graph database contains a pre-computed list of all non-isomorphic
# vertex-transitive graphs for small numbers of vertices. We retrieve this list.
vt_graph_list = graphs.vertex_transitive_graphs_list(N)

# We iterate through the list of graphs and count them based on their degree.
for g in vt_graph_list:
    # A vertex-transitive graph is also a regular graph, so all its vertices
    # have the same degree. We can find the degree of the graph by getting
    # the degree of any single vertex (e.g., the first one).
    # The list includes the empty graph on 8 vertices, which has degree 0.
    degree = g.degree(g.vertices()[0])
    n_j[degree] += 1

# The problem asks for the output in the format [n_0, n_1, ..., n_7].
# The list n_j now contains these counts.
# We format the string explicitly to ensure each number is present.
result_string = f"[{n_j[0]}, {n_j[1]}, {n_j[2]}, {n_j[3]}, {n_j[4]}, {n_j[5]}, {n_j[6]}, {n_j[7]}]"

print("The numbers n_j of isomorphism classes of vertex-transitive graphs with 8 vertices and vertex degree precisely j are:")
print(result_string)

# The final answer is provided in the specified format.
print(f"<<<{result_string}>>>")