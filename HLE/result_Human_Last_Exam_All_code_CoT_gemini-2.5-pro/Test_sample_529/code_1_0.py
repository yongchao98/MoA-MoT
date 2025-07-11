import collections

def solve_tractable_variant(graph, k):
    """
    Solves a tractable variant of the CountAns problem.

    This function calculates the sum of the k-th powers of the degrees of all 
    vertices in the graph. This corresponds to counting the total number of pairs 
    ((v_1, ..., v_k), y) where y is a common neighbor of v_1, ..., v_k.

    Args:
        graph (dict): The graph represented as an adjacency list.
                      Example: {'a': ['b', 'c'], 'b': ['a'], 'c': ['a']}
        k (int): The integer parameter from the problem description.
    """
    if not isinstance(k, int) or k <= 0:
        print("Error: k must be a positive integer.")
        return

    if not graph:
        print("The graph is empty. Result is 0.")
        return
        
    degrees = {node: len(neighbors) for node, neighbors in graph.items()}
    
    # Ensure all nodes mentioned as neighbors are in the degrees dict
    all_nodes = set(graph.keys())
    for neighbors in graph.values():
        for neighbor in neighbors:
            all_nodes.add(neighbor)
    
    for node in all_nodes:
        if node not in degrees:
            degrees[node] = 0

    # To handle disconnected components correctly, we need the degree of every node
    # in the graph, not just those with outgoing edges. Let's recalculate degrees
    # from a list of all edges.
    all_nodes = set(graph.keys())
    edge_list = []
    for node, neighbors in graph.items():
        for neighbor in neighbors:
            all_nodes.add(neighbor)
            edge_list.append(tuple(sorted((node, neighbor))))

    # Recalculate degrees properly for an undirected graph
    final_degrees = collections.defaultdict(int)
    for u, v in set(edge_list):
        final_degrees[u] += 1
        final_degrees[v] += 1
    
    # Ensure nodes with degree 0 are included
    for node in all_nodes:
        if node not in final_degrees:
            final_degrees[node] = 0
            
    total_count = 0
    equation_parts = []

    # Sort nodes for a consistent output order
    sorted_nodes = sorted(final_degrees.keys())

    for node in sorted_nodes:
        degree = final_degrees[node]
        term = degree ** k
        total_count += term
        equation_parts.append(f"{degree}^{k}")

    equation = " + ".join(equation_parts)
    print(f"The equation is: {equation}")
    print(f"The total count is: {total_count}")


# Example Usage:
# A graph with 4 vertices forming a star shape with 'a' at the center.
# 'a' is connected to 'b', 'c', 'd'.
G = {
    'a': ['b', 'c', 'd'],
    'b': ['a'],
    'c': ['a'],
    'd': ['a']
}
k_param = 3

solve_tractable_variant(G, k_param)