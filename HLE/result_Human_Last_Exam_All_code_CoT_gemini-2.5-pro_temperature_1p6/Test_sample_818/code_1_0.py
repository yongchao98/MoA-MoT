import itertools

def count_slices(num_vertices, edges):
    """
    Counts the number of "slices" for a given graph.
    A slice corresponds to a non-empty, proper subset of edges E1 such that
    the subgraph formed by E1 has all vertex degrees equal to 1 or 2.
    The number of slices N(G) is half the count of such subsets, as E1 and (E \ E1)
    form the same partition.
    """
    num_edges = len(edges)
    valid_subsets_count = 0

    # Iterate through all non-empty, proper subsets of edges
    for i in range(1, 2**num_edges - 1):
        subset_edges = []
        for j in range(num_edges):
            if (i >> j) & 1:
                subset_edges.append(edges[j])

        degrees = [0] * num_vertices
        for u, v in subset_edges:
            degrees[u] += 1
            degrees[v] += 1
        
        # The slice condition requires every vertex to be incident to edges from both classes.
        # This is equivalent to the subgraph G1=(V, E1) having degrees d(v) where 1 <= d(v) <= 2.
        is_valid = all(1 <= d <= 2 for d in degrees)
        
        if is_valid:
            valid_subsets_count += 1
    
    # Each slice is a partition {E1, E2}, but our loop counts E1 and E2 separately.
    return valid_subsets_count // 2

def solve_m_n():
    """
    Calculates M(0), M(3), and M(5) and prints the results.
    """
    # Define the graphs to be tested, ordered by number of vertices.
    # Edge lists for all simple cubic graphs with up to 8 vertices.
    graphs_to_test = {
        4: [
            ([(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)], "K4"),
        ],
        6: [
            ([(0,1), (1,2), (2,0), (3,4), (4,5), (5,3), (0,3), (1,4), (2,5)], "Prism"),
            ([(0,3), (0,4), (0,5), (1,3), (1,4), (1,5), (2,3), (2,4), (2,5)], "K3,3"),
        ],
        8: [
            # Cube graph
            ([(0,1),(1,2),(2,3),(3,0), (4,5),(5,6),(6,7),(7,4), (0,4),(1,5),(2,6),(3,7)], "Cube"),
            # The other four 8-vertex cubic graphs
            ([(0,1),(0,2),(0,3),(1,2),(1,4),(2,5),(3,6),(3,7),(4,6),(4,7),(5,6),(5,7)], "G8_2"),
            ([(0,1),(0,2),(0,3),(1,4),(1,5),(2,4),(2,6),(3,5),(3,7),(4,7),(5,6),(6,7)], "G8_3"),
            ([(0,1),(0,2),(0,3),(1,2),(1,4),(2,5),(3,6),(3,7),(4,5),(4,6),(5,7),(6,7)], "G8_4"),
            ([(0,1),(0,5),(0,6),(1,2),(1,7),(2,3),(2,4),(3,5),(3,7),(4,5),(4,6),(6,7)], "G8_5"),
        ],
    }

    # M(0)
    m0 = "none"
    print("For M(0), we need N(G) = 0. It is a known result that for any simple cubic graph G, N(G) > 0. So no such graph exists.")
    print(f"M(0) = {m0}")

    # M(3)
    m3 = "none"
    for m in sorted(graphs_to_test.keys()):
        for edges, name in graphs_to_test[m]:
            n_g = count_slices(m, edges)
            if n_g % 3 == 0:
                m3 = m
                print(f"For M(3), found graph {name} on {m} vertices with N(G) = {n_g}, which is a multiple of 3.")
                break
        if m3 != "none":
            break
    print(f"M(3) = {m3}")
    
    # M(5)
    m5 = "none"
    for m in sorted(graphs_to_test.keys()):
        for edges, name in graphs_to_test[m]:
            n_g = count_slices(m, edges)
            if n_g % 5 == 0 and n_g > 0:
                m5 = m
                print(f"For M(5), found graph {name} on {m} vertices with N(G) = {n_g}, which is a multiple of 5.")
                break
        if m5 != "none":
            break
    print(f"M(5) = {m5}")

    # Final answer in the required format
    final_answer = f"{m0},{m3},{m5}"
    print("Final answer format M(0),M(3),M(5):")
    print(final_answer)
    return final_answer

if __name__ == '__main__':
    final_result = solve_m_n()
    print(f'<<<{final_result}>>>')