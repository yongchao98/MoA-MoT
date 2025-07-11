import itertools

def solve_and_print():
    """
    Calculates the number of specified torsion elements by analyzing
    the E8 Dynkin diagram.
    """
    # Represent the E8 graph. The nodes are numbered 1 to 8.
    # Diagram: 1-2-3-4-5-6-7, with node 8 attached to node 3.
    nodes = set(range(1, 9))
    edges = {(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (3, 8)}

    adj = {i: set() for i in nodes}
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)

    # Find all induced A4 subgraphs (paths on 4 vertices).
    a4_subgraphs = []
    for j_nodes_tuple in itertools.combinations(nodes, 4):
        j_nodes = set(j_nodes_tuple)
        
        degrees_in_subgraph = {v: 0 for v in j_nodes}
        num_edges_in_subgraph = 0
        
        for u, v in itertools.combinations(j_nodes, 2):
            if v in adj[u]:
                degrees_in_subgraph[u] += 1
                degrees_in_subgraph[v] += 1
                num_edges_in_subgraph +=1
        
        # A path graph on 4 vertices has 3 edges and a degree sequence of (1, 2, 2, 1).
        if num_edges_in_subgraph == 3:
            if sorted(degrees_in_subgraph.values()) == [1, 1, 2, 2]:
                a4_subgraphs.append(j_nodes)
    
    # For each A4, count valid disjoint A1 choices.
    total_count = 0
    calculation_steps = []
    
    # Sort for a deterministic and readable output.
    sorted_a4_subgraphs = sorted([tuple(sorted(list(s))) for s in a4_subgraphs])
    
    print("The total number of elements is found by identifying all valid pairs of subgraphs (A4, A1).")
    print("There are 6 distinct subgraphs of type A4 in the E8 diagram.")
    print("For each A4 subgraph, we count the number of available vertices for a disjoint A1 component:\n")
    
    a4_details = []
    for a4_nodes_tuple in sorted_a4_subgraphs:
        a4_nodes = set(a4_nodes_tuple)
        
        # Find neighbors of the A4 subgraph to identify vertices that cannot be used.
        neighbors = set()
        for u in a4_nodes:
            neighbors.update(adj[u])
        
        forbidden = a4_nodes.union(neighbors)
        a1_choices = sorted(list(nodes.difference(forbidden)))
        
        count = len(a1_choices)
        a4_details.append(f"For the A4 on vertices {a4_nodes}, there are {count} choices for A1: {a1_choices}")
        calculation_steps.append(count)
        total_count += count

    for detail in a4_details:
        print(detail)

    # Print the final calculation as a sum.
    equation_str = " + ".join(map(str, calculation_steps))
    print(f"\nThe total number of such elements is the sum of these counts:")
    print(f"Total = {equation_str} = {total_count}")

solve_and_print()