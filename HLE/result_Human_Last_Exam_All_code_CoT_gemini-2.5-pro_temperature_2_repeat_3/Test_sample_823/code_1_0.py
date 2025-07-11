import networkx as nx

def demonstrate_counterexample_for_E():
    """
    This script demonstrates that statement E is false using a counterexample.
    Statement E: For a class of graphs C with bounded degree and unbounded 
                 treewidth, the graphs in C contain clique-minors of unbounded size.
                 
    We will use the class of n x n grid graphs as our C.
    This class satisfies the premises:
    1. Bounded degree: The maximum degree is 4 for any n >= 2.
    2. Unbounded treewidth: The treewidth of an n x n grid is n.
    
    We will show that this class does not satisfy the conclusion of statement E.
    """
    
    print("--- Demonstrating a counterexample for statement E ---")
    print("Let C be the class of n x n grid graphs.")
    print("This class has bounded degree (<= 4) and unbounded treewidth (tw(G_n) = n).")
    print("We check if this class must have clique-minors of unbounded size.\n")
    
    # Iterate through grids of increasing size
    for n in range(2, 7):
        G = nx.grid_2d_graph(n, n)
        
        # Treewidth of an n x n grid is n. This is a known result.
        treewidth = n
        
        # Grid graphs are planar. We can verify this property.
        is_planar, _ = nx.check_planarity(G)
        
        print(f"Analysis of the {n}x{n} grid:")
        
        # Print the relevant numbers and properties for the grid.
        # The equation here is an informal one: Treewidth(G_n) = n
        print(f"  Treewidth = {treewidth}")
        print(f"  Planar = {is_planar}")
        
        # According to Wagner's Theorem, a graph is planar if and only if
        # it does not have K_5 (clique of size 5) or K_{3,3} as a minor.
        # This means the largest clique minor any planar graph can have is K_4.
        max_clique_minor_size = 4
        print(f"  The size of the largest clique-minor is bounded by {max_clique_minor_size}.")
        print("")

    print("--- Conclusion ---")
    print("The class of grid graphs has treewidth that grows with n, so it's unbounded.")
    print("However, the size of the largest clique-minor is always bounded by 4.")
    print("This shows that statement E is false.")

demonstrate_counterexample_for_E()