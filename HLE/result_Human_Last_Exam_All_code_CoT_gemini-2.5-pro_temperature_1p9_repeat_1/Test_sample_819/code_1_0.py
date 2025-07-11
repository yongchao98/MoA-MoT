import networkx as nx

def find_simple_paths_of_length_k(graph, k):
    """
    Finds all simple paths of a specific length k in a graph.
    A simple path is one with no repeated vertices.

    Args:
        graph (nx.Graph): The graph to search in.
        k (int): The desired length of the paths (number of edges).

    Returns:
        A list of paths, where each path is a list of vertices.
    """
    paths = []
    for source_node in graph.nodes():
        # Target nodes must be different from source_node
        target_nodes = [n for n in graph.nodes() if n != source_node]
        # networkx's all_simple_paths finds paths of any length. We filter them.
        for path in nx.all_simple_paths(graph, source=source_node, target=target_nodes, cutoff=k):
            if len(path) == k + 1:
                paths.append(path)
    return paths

def count_query_answers(graph):
    """
    Counts the number of answers for the query phi(x,y) asking for a
    simple path of length 3 between x and y.

    An answer is an ordered pair (x, y) of vertices. We count the number
    of unique pairs.
    """
    path_length = 3
    simple_paths = find_simple_paths_of_length_k(graph, path_length)

    # An answer is a tuple (start_node, end_node)
    # We use a set to count unique answers
    answers = set()
    for path in simple_paths:
        start_node = path[0]
        end_node = path[-1]
        answers.add((start_node, end_node))

    return len(answers), sorted(list(answers))

# Create graph G1: two disjoint triangles (2C3)
g1 = nx.Graph()
g1.add_edges_from([(1, 2), (2, 3), (3, 1)]) # First C3
g1.add_edges_from([(4, 5), (5, 6), (6, 4)]) # Second C3

# Create graph G2: a 6-cycle (C6)
g2 = nx.cycle_graph(6)

# Calculate and print the number of answers for both graphs
num_answers_g1, answers_g1 = count_query_answers(g1)
print(f"Graph G1 (2*C3):")
print(f"The number of answers is: {num_answers_g1}")
print(f"The answers (pairs (x,y)) are: {answers_g1}")

print("-" * 20)

num_answers_g2, answers_g2 = count_query_answers(g2)
print(f"Graph G2 (C6):")
print(f"The number of answers is: {num_answers_g2}")
print(f"The answers (pairs (x,y)) are: {answers_g2}")

if num_answers_g1 != num_answers_g2:
    print("\nConclusion: The number of answers is different for the two graphs.")
else:
    print("\nConclusion: The number of answers is the same for the two graphs.")