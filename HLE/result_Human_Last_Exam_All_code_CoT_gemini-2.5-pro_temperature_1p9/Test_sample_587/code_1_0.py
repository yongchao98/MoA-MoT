import collections

def check_degree_at_most_2(graph_h):
    """
    Checks if all vertices in a graph have a degree of at most 2.

    The graph is represented as an adjacency list (dictionary of lists).
    A graph where all vertices have degree at most 2 is a collection of
    disjoint paths and cycles. Such graphs have a treewidth of at most 2.
    """
    if not graph_h:
        return True # An empty graph has degrees at most 2.
        
    for vertex in graph_h:
        degree = len(graph_h[vertex])
        print(f"Vertex {vertex} has degree {degree}.")
        if degree > 2:
            print(f"Vertex {vertex} has degree greater than 2.")
            return False
    print("All vertices have a degree of at most 2.")
    return True

def main():
    """
    Main function to demonstrate the degree check on sample graphs.
    """
    # H1: A 4-cycle (C4). Degrees are all 2.
    # This represents a class H of graphs with degree at most 2.
    # According to statement C, counting such subgraphs could be #W[1]-complete.
    # However, since its treewidth is bounded, the problem is FPT, making C false.
    h1 = {
        'a': ['b', 'd'],
        'b': ['a', 'c'],
        'c': ['b', 'd'],
        'd': ['a', 'c']
    }
    print("--- Analyzing graph H1 (a 4-cycle) ---")
    result1 = check_degree_at_most_2(h1)
    print(f"Does H1 have maximum degree at most 2? {result1}\n")

    # H2: A 4-clique (K4). Degrees are all 3.
    # This graph does not have degree at most 2. It has unbounded treewidth
    # (as part of the class of cliques), which leads to hardness.
    h2 = {
        1: [2, 3, 4],
        2: [1, 3, 4],
        3: [1, 2, 4],
        4: [1, 2, 3]
    }
    print("--- Analyzing graph H2 (a 4-clique) ---")
    result2 = check_degree_at_most_2(h2)
    print(f"Does H2 have maximum degree at most 2? {result2}\n")

    print("Explanation:")
    print("Statement C suggests there's a class H of graphs with degree at most 2 for which subgraph counting is #W[1]-complete.")
    print("As the code illustrates, graphs with degree at most 2 (like H1) are simple collections of paths and cycles.")
    print("These graphs have a bounded treewidth (at most 2).")
    print("It is a known theorem that for any class H with bounded treewidth, the subgraph counting problem is Fixed-Parameter Tractable (FPT).")
    print("Therefore, statement C is false.")


if __name__ == "__main__":
    main()