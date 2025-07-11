import itertools

def solve_continuum_decomposition():
    """
    Solves a graph theory problem equivalent to the topological puzzle.

    The problem is to find the maximum number of edges in a graph on 5 vertices
    such that it does not contain a K_3 (triangle) and a disjoint K_2 (edge)
    as a subgraph.
    """
    vertices = list(range(5))
    all_edges = list(itertools.combinations(vertices, 2))

    # Generate all 10 forbidden sets of edges (K_3 U K_2)
    forbidden_sets = []
    for tri_verts in itertools.combinations(vertices, 3):
        remaining_verts = tuple(v for v in vertices if v not in tri_verts)
        
        # The three edges of the triangle
        tri_edges = set(itertools.combinations(tri_verts, 2))
        
        # The one disjoint edge
        disjoint_edge = (remaining_verts[0], remaining_verts[1])

        # A forbidden set has 4 edges
        forbidden_set = tri_edges.copy()
        forbidden_set.add(disjoint_edge)
        forbidden_sets.append(forbidden_set)

    # Search for the largest valid subset of edges, from k=10 down to 1
    for k in range(len(all_edges), 0, -1):
        for edge_subset in itertools.combinations(all_edges, k):
            current_set = set(edge_subset)
            is_valid = True
            for f_set in forbidden_sets:
                if f_set.issubset(current_set):
                    is_valid = False
                    break
            
            if is_valid:
                print(f"The largest number n is {k}.")
                print("\nThis value is found by identifying the problem as finding the maximum number of edges in a 5-vertex graph that does not contain a triangle and a disjoint edge. The code confirmed this number by checking all possibilities.")
                
                print("\nThe decomposition equation with this n is:")
                # As requested, printing each number in the final equation.
                equation_parts = [f"A_{i+1}" for i in range(k)]
                final_equation = "X = " + " U ".join(equation_parts)
                print(final_equation)
                return

solve_continuum_decomposition()
