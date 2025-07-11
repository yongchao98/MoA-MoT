import networkx as nx
from itertools import product

def construct_rooks_graph():
    """Constructs the L(K_4,4) graph, also known as the 4x4 Rook's graph."""
    G = nx.Graph()
    nodes = list(product(range(4), range(4)))
    G.add_nodes_from(nodes)
    for u in nodes:
        for v in nodes:
            if u == v:
                continue
            # Two vertices are adjacent if they share a row or a column
            if u[0] == v[0] or u[1] == v[1]:
                G.add_edge(u, v)
    return G

def construct_shrikhande_graph():
    """Constructs the Shrikhande graph as a Cayley graph on Z_4 x Z_4."""
    G = nx.Graph()
    nodes = list(product(range(4), range(4)))
    G.add_nodes_from(nodes)
    # Connection set for the Cayley graph
    connection_set = {(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)}
    for u in nodes:
        for v in nodes:
            if u == v:
                continue
            # Check if the difference (mod 4) is in the connection set
            diff = ((u[0] - v[0]) % 4, (u[1] - v[1]) % 4)
            if diff in connection_set:
                G.add_edge(u, v)
    return G

def count_five_cycles(graph):
    """Counts the number of simple cycles of length 5 in a graph."""
    # nx.simple_cycles finds all simple cycles (as lists of nodes).
    # We filter for those with length 5 and sum them up.
    return sum(1 for cycle in nx.simple_cycles(graph) if len(cycle) == 5)

def main():
    """
    Main function to construct graphs, count cycles, and print the results.
    """
    print("This script checks if two Strongly Regular Graphs with the same parameters")
    print("can have a different number of 5-cycles.\n")
    
    # The parameters for both graphs are (n, d, lambda, mu)
    n = 16
    d = 6
    l = 2
    m = 2
    
    print(f"The chosen SRG parameters are (n, d, lambda, mu) = ({n}, {d}, {l}, {m}).\n")
    
    # Construct the two graphs
    rooks_graph = construct_rooks_graph()
    shrikhande_graph = construct_shrikhande_graph()
    
    # Count the 5-cycles in each graph
    rooks_cycles = count_five_cycles(rooks_graph)
    shrikhande_cycles = count_five_cycles(shrikhande_graph)
    
    # Output the results
    print("Graph 1: Rook's Graph L(K_4,4)")
    print(f"Number of 5-cycles = {rooks_cycles}")
    
    print("\nGraph 2: Shrikhande Graph")
    print(f"Number of 5-cycles = {shrikhande_cycles}\n")
    
    if rooks_cycles != shrikhande_cycles:
        print("Conclusion: Yes, two SRGs with the same parameters can have a different number of 5-cycles.")
    else:
        print("Conclusion: In this case, the number of 5-cycles is the same.")

if __name__ == "__main__":
    main()
<<<Yes>>>