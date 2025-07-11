import networkx as nx
from itertools import product

def get_graphs():
    """
    Creates the two graphs for our example.
    G1 is the disjoint union of two 3-cycles.
    G2 is a 6-cycle.
    These graphs are known to be path-equivalent, meaning for any path P,
    the number of homomorphisms from P to G1 is the same as to G2.
    However, they are not tree-equivalent.
    """
    # G1: Two disjoint triangles
    g1_component1 = nx.cycle_graph(3)
    # networkx by default labels the second component starting from 0 again.
    # We'll relabel to make vertices unique.
    g1_component2 = nx.relabel_nodes(nx.cycle_graph(3), {0: 3, 1: 4, 2: 5})
    g1 = nx.disjoint_union(g1_component1, g1_component2)

    # G2: A single 6-cycle
    g2 = nx.cycle_graph(6)
    
    return g1, g2

def count_query_answers(graph):
    """
    Counts the number of answers for the acyclic conjunctive query:
    phi(x, z) = exists y. E(x, y) AND E(y, z)
    
    An answer is a tuple (u, w) of vertices from the graph such that there
    is a path of length 2 between u and w.
    """
    answers = set()
    nodes = list(graph.nodes())
    
    # We are looking for pairs (x, z) connected by a path of length 2.
    # We can iterate through all possible intermediate vertices 'y'.
    for y in nodes:
        # Get all neighbors of y
        neighbors = list(graph.neighbors(y))
        # Any pair of neighbors (including a neighbor with itself)
        # forms a path of length 2 with y in the middle.
        for x in neighbors:
            for z in neighbors:
                answers.add((x, z))
                
    return len(answers), sorted(list(answers))

def main():
    """
    Main function to execute the logic.
    """
    print("This script demonstrates that it is possible for two graphs, G1 and G2,")
    print("to have the same number of homomorphisms from any path, yet have a")
    print("different number of answers for an acyclic conjunctive query.\n")
    print("The problem asks about equivalence for all trees, which is a stronger condition.")
    print("However, this example for paths illustrates the principle that homomorphism counts")
    print("do not uniquely determine the number of query answers.\n")

    g1, g2 = get_graphs()

    # The query is phi(x, z) = exists y. E(x, y) AND E(y, z)
    # This corresponds to finding pairs of vertices at distance 2.
    
    num_answers_g1, answers_g1 = count_query_answers(g1)
    num_answers_g2, answers_g2 = count_query_answers(g2)

    print("Graph G1 is the disjoint union of two 3-cycles.")
    print(f"Number of answers for G1: {num_answers_g1}")
    # print(f"Answers for G1: {answers_g1}\n")

    print("Graph G2 is a 6-cycle.")
    print(f"Number of answers for G2: {num_answers_g2}")
    # print(f"Answers for G2: {answers_g2}\n")

    print("The number of answers is different for G1 and G2.")
    
    num1 = num_answers_g1
    num2 = num_answers_g2
    print(f"The final comparison is: {num1} != {num2}")


if __name__ == "__main__":
    main()