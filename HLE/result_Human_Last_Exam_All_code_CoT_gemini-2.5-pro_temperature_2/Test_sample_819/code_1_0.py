import itertools

def find_answers(graph_edges, graph_nodes):
    """
    Finds the number of answers for the acyclic conjunctive query:
    phi(x) = Exists u,v,w. E(x,u) AND E(u,v) AND E(v,w) AND E(w,u)
    
    An answer is a vertex x that is a neighbor of a vertex u,
    where u is part of a triangle {u,v,w}.
    """
    answers = set()
    
    # Create an adjacency set for faster lookups
    adj = {node: set() for node in graph_nodes}
    for n1, n2 in graph_edges:
        adj[n1].add(n2)
        adj[n2].add(n1)

    # Iterate through all combinations of 3 vertices to find triangles
    for u, v, w in itertools.permutations(graph_nodes, 3):
        # Check if {u, v, w} forms a triangle
        if v in adj[u] and w in adj[v] and u in adj[w]:
            # u,v,w form a triangle.
            # Now, find all neighbors of u. These are potential answers 'x'.
            for x in adj[u]:
                answers.add(x)

    return len(answers)

def solve():
    """
    Solves the problem by constructing the graphs G1 and G2,
    and then finding the number of answers to the query for each.
    """
    # G1: Cycle graph on 6 vertices (C6)
    g1_nodes = [0, 1, 2, 3, 4, 5]
    g1_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
    
    # G2: Disjoint union of two triangles (2*K3)
    g2_nodes = [0, 1, 2, 3, 4, 5]
    g2_edges = [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)]

    # The query graph for phi(x) is acyclic (a "paw" graph).
    # Its variables are x, u, v, w.
    # Its atoms (edges) are (x,u), (u,v), (v,w), (w,u).
    # The vertices v,w,u form a triangle, and x is a neighbor of u.

    num_answers_g1 = find_answers(g1_edges, g1_nodes)
    num_answers_g2 = find_answers(g2_edges, g2_nodes)
    
    print("G1 is a cycle on 6 vertices (C6).")
    print("G2 is the disjoint union of two triangles (2*K3).")
    print("For any tree T, the number of homomorphisms from T to G1 is equal to the number from T to G2.\n")
    print("Let the acyclic query be phi(x) = Exists u,v,w such that (x,u), (u,v), (v,w), (w,u) are all edges.")
    print("An answer is a vertex 'x' that is adjacent to a vertex 'u' which is part of a triangle.\n")

    print(f"In G1 (C6), a triangle does not exist, so there are no vertices 'u', 'v', 'w' satisfying the condition.")
    print(f"Number of answers for phi in G1: {num_answers_g1}")

    print(f"\nIn G2 (2*K3), triangles exist. For any vertex 'u' in a triangle, its neighbors 'x' within that triangle are answers.")
    print(f"Number of answers for phi in G2: {num_answers_g2}")

    if num_answers_g1 != num_answers_g2:
        print("\nThe number of answers is different for the two graphs.")
    else:
        print("\nThe number of answers is the same for the two graphs.")

solve()