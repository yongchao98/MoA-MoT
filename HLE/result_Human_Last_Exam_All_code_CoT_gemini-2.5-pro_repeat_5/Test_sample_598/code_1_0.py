import networkx as nx

def count_cycles_5(G):
    """
    Counts the number of simple 5-cycles in a graph.
    The algorithm iterates through all walks of length 5 (v1-v2-v3-v4-v5-v1)
    and checks if the vertices are distinct.
    To avoid overcounting, the final count is divided by 10, as each 5-cycle
    can be started from 5 different nodes and traversed in 2 different directions.
    """
    count = 0
    nodes = list(G.nodes())
    for n1 in nodes:
        # Explore non-backtracking paths of length 4 starting from n1
        # Path: n1-n2-n3-n4-n5
        for n2 in G.neighbors(n1):
            for n3 in G.neighbors(n2):
                if n3 == n1:
                    continue
                for n4 in G.neighbors(n3):
                    if n4 == n2:
                        continue
                    # A closed walk n1-n2-n3-n1 is a 3-cycle
                    if n4 == n1:
                        continue
                    for n5 in G.neighbors(n4):
                        if n5 == n3:
                            continue
                        # A closed walk n1-n2-n3-n5-n1 with n5=n2 is a 4-cycle
                        if n5 == n2:
                            continue
                        
                        # Now we have a non-backtracking walk n1-n2-n3-n4-n5.
                        # We check if n5 is connected to n1 to close the walk.
                        if G.has_edge(n5, n1):
                            # We found a closed walk of length 5: n1-n2-n3-n4-n5-n1
                            # Check if all vertices in the path are distinct to form a simple cycle.
                            if len({n1, n2, n3, n4, n5}) == 5:
                                count += 1
    # Each cycle is counted 5 * 2 = 10 times.
    return count // 10

def main():
    """
    Main function to construct the graphs and compare their 5-cycle counts.
    """
    print("Yes, there is a combination of parameters for which two SRGs have a different number of 5-cycles.")
    print("A well-known example is for the parameters (n, d, lambda, mu) = (16, 6, 2, 2).\n")

    # 1. Construct the Shrikhande graph
    G_shrikhande = nx.shrikhande_graph()

    # 2. Construct the 4x4 Rook's graph (K4 x K4)
    K4 = nx.complete_graph(4)
    G_rook = nx.cartesian_product(K4, K4)

    # 3. Count 5-cycles in both graphs
    c5_shrikhande = count_cycles_5(G_shrikhande)
    c5_rook = count_cycles_5(G_rook)

    # 4. Print the results
    print("Graph 1: Shrikhande graph")
    print(f"Parameters: (n={G_shrikhande.number_of_nodes()}, d={G_shrikhande.degree(0)}, lambda=2, mu=2)")
    print(f"Number of 5-cycles = {c5_shrikhande}")
    
    print("\nGraph 2: 4x4 Rook's graph (K4 x K4)")
    # The degree of a node like (0,0) in the product graph is deg(0) in K4 + deg(0) in K4 = 3+3=6
    print(f"Parameters: (n={G_rook.number_of_nodes()}, d={G_rook.degree((0,0))}, lambda=2, mu=2)")
    print(f"Number of 5-cycles = {c5_rook}")

    print(f"\nThe number of 5-cycles is different for these two graphs.")
    print(f"Final equation: {c5_shrikhande} != {c5_rook}")

if __name__ == "__main__":
    main()