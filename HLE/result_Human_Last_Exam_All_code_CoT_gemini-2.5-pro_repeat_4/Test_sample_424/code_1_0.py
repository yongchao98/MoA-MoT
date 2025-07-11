import networkx as nx

def solve_planar_set_problem():
    """
    Solves the problem by modeling the planar set as a graph and analyzing its cut points.
    """
    # 1. Define the graph nodes (junctions and endpoints) and edges (segments and arcs)
    nodes = [
        'p(0, 1/2)', 'p(0, 3/2)', 'p(1/2, 0)', 'p(-3/2, 0)',
        'p(-1/2, 0)', 'p(0, -1/2)', 'p(-1/2, 1)', 'p(1/2, 1)',
        'j(0, 1)', 'j(1, 0)', 'j(-1, 0)', 'j(0, -1)',
        'j(3/2, 0)', 'j(0, -3/2)'
    ]

    edges = [
        # Unit circle arcs between junction points
        ('j(0, 1)', 'j(1, 0)'), ('j(1, 0)', 'j(0, -1)'), 
        ('j(0, -1)', 'j(-1, 0)'), ('j(-1, 0)', 'j(0, 1)'),
        # Line segment {0} x [1/2, 3/2]
        ('p(0, 1/2)', 'j(0, 1)'), ('j(0, 1)', 'p(0, 3/2)'),
        # Line segment [1/2, 3/2] x {0}
        ('p(1/2, 0)', 'j(1, 0)'), ('j(1, 0)', 'j(3/2, 0)'),
        # Line segment [-3/2, -1/2] x {0}
        ('p(-3/2, 0)', 'j(-1, 0)'), ('j(-1, 0)', 'p(-1/2, 0)'),
        # Line segment {0} x [-3/2, -1/2]
        ('p(0, -1/2)', 'j(0, -1)'), ('j(0, -1)', 'j(0, -3/2)'),
        # Line segment [-1/2, 1/2] x {1}
        ('p(-1/2, 1)', 'j(0, 1)'), ('j(0, 1)', 'p(1/2, 1)'),
        # Bottom-right quarter circle of radius 3/2
        ('j(3/2, 0)', 'j(0, -3/2)')
    ]

    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    # 2. Identify junction points (degree > 2) and leaf nodes (degree == 1)
    junction_points = [node for node, degree in G.degree() if degree > 2]
    leaf_nodes = {node for node, degree in G.degree() if degree == 1}

    solution_points = []
    
    print("Analyzing junction points to find which ones are cut points that create >= 3 components:")
    # 3. Analyze each junction point
    for j_node in junction_points:
        neighbors = list(G.neighbors(j_node))
        
        # 4. Count dangling branches connected to the junction
        dangling_branches = 0
        for neighbor in neighbors:
            if neighbor in leaf_nodes:
                dangling_branches += 1
        
        # The number of components is 1 (the main graph body) + the number of separated dangling branches
        num_components = dangling_branches + 1
        
        print(f"\n- Removing point {j_node}:")
        print(f"  - It connects to {len(neighbors)} branches in the graph model.")
        print(f"  - Of these, {dangling_branches} are dangling branches (leading to a free end).")
        print(f"  - The number of resulting components = (dangling branches) + 1 = {dangling_branches} + 1 = {num_components}")

        if num_components >= 3:
            solution_points.append(j_node)
            print(f"  - Verdict: This point's removal creates {num_components} (>= 3) components. It is a solution.")
        else:
            print(f"  - Verdict: This point's removal creates {num_components} (< 3) components. Not a solution.")

    print("\n" + "="*50)
    print("Final Conclusion:")
    print(f"The points whose removal results in 3 or more components are: {solution_points}")
    print(f"The total number of such points is: {len(solution_points)}")

solve_planar_set_problem()