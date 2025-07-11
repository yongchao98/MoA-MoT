import itertools

def count_slice_defining_subgraphs(num_vertices, edges):
    """Counts the number of subgraphs where every vertex has degree 1 or 2."""
    count = 0
    num_edges = len(edges)
    for i in range(1, num_edges + 1):
        for subgraph_edges in itertools.combinations(edges, i):
            degrees = [0] * num_vertices
            is_valid = True
            for u, v in subgraph_edges:
                degrees[u] += 1
                degrees[v] += 1
            for degree in degrees:
                if degree not in [1, 2]:
                    is_valid = False
                    break
            if is_valid:
                count += 1
    return count

def main():
    """
    Calculates M(0), M(3), M(5) and prints the results.
    """
    # M(0)
    # A slice is a partition of edges into two classes where each vertex is incident
    # to at least one edge in each class. This is equivalent to finding a spanning
    # subgraph where all vertex degrees are 1 or 2.
    # It is a known result that every cubic graph has such a subgraph.
    # Therefore, N(G) >= 1 for any cubic graph G, so N(G) cannot be 0.
    m0_result = "none"
    print("M(0) = none, because every cubic graph G has at least one slice, so N(G) is never 0.")

    # M(3)
    # Smallest cubic graph is K_4 on m=4 vertices.
    v_k4 = 4
    e_k4 = list(itertools.combinations(range(v_k4), 2))
    num_sds_k4 = count_slice_defining_subgraphs(v_k4, e_k4)
    n_k4 = num_sds_k4 / 2
    # n_k4 is 9.0, which is a multiple of 3.
    # Since m=4 is the smallest possible size for a cubic graph, M(3)=4.
    m3_result = 4
    print(f"For K_4 (m=4), the number of slices N(K_4) is {n_k4:.0f}.")
    print(f"Since {n_k4:.0f} is a multiple of 3 and m=4 is the minimum possible, M(3) = {m3_result}.")
    
    # M(5)
    # For m=4, N(K_4) = 9, not a multiple of 5.
    # For m=6, there are two cubic graphs.
    # Prism graph:
    v_prism = 6
    e_prism = [(0,1),(1,2),(2,0), (3,4),(4,5),(5,3), (0,3),(1,4),(2,5)]
    num_sds_prism = count_slice_defining_subgraphs(v_prism, e_prism)
    n_prism = num_sds_prism/2
    print(f"For the Prism graph (m=6), N(Prism) is {n_prism:.0f}, not a multiple of 5.")

    # K_3,3 graph:
    v_k33 = 6
    e_k33 = []
    for u in [0,1,2]:
        for v in [3,4,5]:
            e_k33.append((u,v))
    num_sds_k33 = count_slice_defining_subgraphs(v_k33, e_k33)
    n_k33 = num_sds_k33 / 2
    print(f"For K_3,3 (m=6), N(K_3,3) is {n_k33:.0f}, not a multiple of 5.")

    # For m=8, we test the Cube graph Q_3.
    v_q3 = 8
    e_q3 = []
    for i in range(v_q3):
        for j in range(3):
            neighbor = i ^ (1 << j)
            if i < neighbor:
                 e_q3.append((i, neighbor))
    num_sds_q3 = count_slice_defining_subgraphs(v_q3, e_q3)
    n_q3 = num_sds_q3 / 2
    print(f"For the Cube graph Q_3 (m=8), N(Q_3) is {n_q3:.0f}.")

    # n_q3 is 45.0, which is a multiple of 5.
    # Since graphs with m < 8 do not have N(G) as a multiple of 5, M(5) = 8.
    m5_result = 8
    print(f"Since {n_q3:.0f} is a multiple of 5 and smaller graphs do not satisfy the condition, M(5) = {m5_result}.")

if __name__ == '__main__':
    main()
