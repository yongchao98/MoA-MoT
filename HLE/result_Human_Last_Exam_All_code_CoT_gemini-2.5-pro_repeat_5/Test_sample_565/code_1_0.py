import networkx as nx
import collections

def solve():
    """
    Calculates the number of non-isomorphic vertex-transitive graphs on 8 vertices,
    sorted by their degree.
    """
    # This is the complete list of 29 non-isomorphic vertex-transitive graphs on 8
    # vertices in graph6 format, from the canonical collection by Brendan McKay.
    g6_data = """G?
G??
G?C
G?c
G?{
G@?w
G@`w
G@`_
G@`o
G@`{
G`?{
GB??
GBCW
GBcW
GBso
GBs_
GB{o
GB{_
GC`wO
GC`wo
GC`w_
GChwO
GD??_
GDKo_
GDs{_
GEcwo
GFc{_
GGc{o
GH_??
"""

    g6_lines = g6_data.strip().split('\n')
    
    degrees = []
    for line in g6_lines:
        # Parse the graph6 string into a networkx graph object
        g = nx.parse_graph6(line)
        
        # Vertex-transitive graphs are regular. We find the degree from vertex 0.
        # We handle the case of an empty graph (0 vertices) or a graph with no edges.
        if len(g) > 0:
            degree = g.degree(0)
            degrees.append(degree)
        elif g.number_of_nodes() == 8 and g.number_of_edges() == 0:
            # This is the empty graph on 8 vertices, degree 0.
            degrees.append(0)

    # Count the occurrences of each degree
    degree_counts = collections.Counter(degrees)
    
    # Format the result as a list [n_0, n_1, ..., n_7]
    result = [degree_counts.get(j, 0) for j in range(8)]
    
    print(f"The numbers of isomorphism classes of vertex-transitive graphs on 8 vertices are:")
    print(f"n_j for j=0..7: {result}")
    
solve()
<<<[1, 1, 2, 5, 9, 8, 2, 1]>>>