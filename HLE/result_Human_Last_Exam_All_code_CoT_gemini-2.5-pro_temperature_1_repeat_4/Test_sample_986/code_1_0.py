import networkx as nx
import itertools

def compute_clique_number_for_finite_case(n):
    """
    This function builds the graph X for a finite poset D_n = {0, ..., n-1}
    and computes its clique number to verify the theoretical result.
    
    Args:
        n (int): The size of the finite set to use.
        
    Returns:
        int: The clique number of the constructed graph.
    """
    # For n<2, the graph is empty. For n=2, it has one vertex.
    # The interesting case is n>=3.
    if n < 3:
        return n-1 if n > 0 else 0

    # The vertices of the graph X are pairs (a, b) from D_n with a < b.
    vertices_X = list(itertools.combinations(range(n), 2))

    # Create the underlying undirected graph of X, let's call it X_u.
    X_u = nx.Graph()
    X_u.add_nodes_from(vertices_X)

    # Add edges to X_u. Two vertices v1=(a,b) and v2=(c,d) are adjacent
    # if the head of one matches the tail of the other.
    for v1, v2 in itertools.combinations(vertices_X, 2):
        a, b = v1
        c, d = v2
        
        # Adjacency is defined by b==c or d==a.
        # The ordering constraints (e.g., a<b<d) are implicitly satisfied
        # by the way the vertices themselves are defined.
        if b == c or d == a:
            X_u.add_edge(v1, v2)

    # Use the networkx library to find the size of the largest clique.
    return nx.graph_clique_number(X_u)

# The theoretical analysis proves the clique number is 2.
# We run a computation on a finite example (n=8) to confirm.
# The result is 2 for any n >= 3.
n_example = 8
clique_number = compute_clique_number_for_finite_case(n_example)

# The final equation is "Clique Number = 2".
# We print the value from our computation.
print("The computed clique number for the example case is:")
print(clique_number)