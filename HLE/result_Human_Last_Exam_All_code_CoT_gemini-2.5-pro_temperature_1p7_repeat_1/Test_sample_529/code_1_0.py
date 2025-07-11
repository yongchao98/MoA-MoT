import networkx as nx
from itertools import product

def count_answers_phi_k(G, k):
    """
    Counts the number of answers for the formula phi_k in graph G.

    An answer is a k-tuple of vertices (x_1, ..., x_k) for which there
    exists a common neighbor y. This algorithm correctly finds the number
    of such unique tuples.

    The time complexity of this algorithm is not fixed-parameter tractable (FPT),
    as the size of the set of answers can be large (up to n^k), and iterating
    through products of neighborhoods is exponential in k.

    Args:
        G (networkx.Graph): The input graph.
        k (int): The parameter k.

    Returns:
        int: The number of answers.
    """
    if not isinstance(k, int) or k < 0:
        raise ValueError("k must be a non-negative integer.")
    if k == 0:
        # The only 0-tuple is the empty tuple, which is satisfied by any graph with at least one vertex
        # as there is no condition to check. Some conventions might return 0 or 1.
        # Assuming we are looking for non-empty tuples if not specified.
        # However, a simple interpretation is that it is satisfied vacuously. Let's return 1 for a non-empty graph.
        return 1 if G.nodes else 0
    if not G.nodes:
        return 0

    answer_tuples = set()
    nodes = list(G.nodes())
    
    # Iterate through all vertices `y` that could be a common neighbor.
    for y in nodes:
        # The neighbors of y are candidates for the answer tuple.
        neighbors = list(G.neighbors(y))
        
        # If y has fewer than 1 neighbor, it cannot be a common neighbor
        # for a k>=1 tuple whose components are from N(y).
        if not neighbors:
            continue
        
        # Any k-tuple formed by the neighbors of y is an answer tuple.
        # itertools.product generates the Cartesian product of the neighbors k times.
        for answer_tuple in product(neighbors, repeat=k):
            answer_tuples.add(answer_tuple)
            
    return len(answer_tuples)

def main():
    # Example Graph
    G = nx.Graph()
    # A 4-star centered at 'c'
    G.add_edges_from([('c', 'L1'), ('c', 'L2'), ('c', 'L3'), ('c', 'L4')])
    # Another node 'x' connected to two leaves of the star
    G.add_edges_from([('x', 'L1'), ('x', 'L2')])
    # An isolated edge
    G.add_edge('a', 'b')

    # Parameter k
    k = 2

    # The number of answers is the size of the union of N(y)^k for all y in V.
    # N(c)^2 -> |{L1,L2,L3,L4}|^2 = 16 answers
    # N(x)^2 -> |{L1,L2}|^2 = 4 answers. These are a subset of the answers from N(c).
    # N(L1)^2 -> |{c,x}|^2 = 4 answers
    # N(L2)^2 -> |{c,x}|^2 = 4 answers. These are the same as from N(L1).
    # N(L3)^2 -> |{c}|^2 = 1 answer, which is ('c', 'c'). This is inside N(L1)^2.
    # N(L4)^2 -> |{c}|^2 = 1 answer, ('c', 'c').
    # N(a)^2 -> |{b}|^2 = 1 answer, ('b', 'b').
    # N(b)^2 -> |{a}|^2 = 1 answer, ('a', 'a').
    
    # The set of all unique answers is the union of:
    # { (Li, Lj) | i,j in {1,2,3,4} }         (16 tuples from N(c))
    # { (v1, v2) | v1,v2 in {c,x} }           (4 tuples from N(L1), N(L2))
    # { ('b','b') }                               (1 tuple from N(a))
    # { ('a','a') }                               (1 tuple from N(b))
    # Total unique answers = 16 + 4 + 1 + 1 = 22.
    
    num_answers = count_answers_phi_k(G, k)

    print(f"Graph details: {len(G.nodes())} nodes, {len(G.edges())} edges.")
    print(f"For parameter k = {k}:")
    print(f"The number of answers of phi_k is: {num_answers}")

if __name__ == "__main__":
    main()