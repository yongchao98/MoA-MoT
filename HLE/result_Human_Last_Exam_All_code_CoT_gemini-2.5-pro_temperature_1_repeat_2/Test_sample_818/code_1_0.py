import itertools

def get_subgraph_degrees(num_vertices, edge_subset):
    """Calculates the degree of each vertex in a subgraph defined by a subset of edges."""
    degrees = {i: 0 for i in range(num_vertices)}
    for u, v in edge_subset:
        degrees[u] += 1
        degrees[v] += 1
    return degrees

def is_slice_defining(num_vertices, edge_subset, total_degrees):
    """Checks if an edge subset defines a slice."""
    if not edge_subset:
        return False
    
    degrees = get_subgraph_degrees(num_vertices, edge_subset)
    
    for v in range(num_vertices):
        # The degree in the subgraph cannot be 0 or the total degree (3 for a cubic graph).
        if degrees[v] == 0 or degrees[v] == total_degrees[v]:
            return False
            
    return True

def count_slices(num_vertices, edges):
    """Counts the number of slices for a given graph."""
    num_slice_defining_sets = 0
    total_degrees = get_subgraph_degrees(num_vertices, edges)

    # Iterate through all non-empty subsets of edges
    for i in range(1, 1 << len(edges)):
        subset = []
        for j in range(len(edges)):
            if (i >> j) & 1:
                subset.append(edges[j])
        
        if is_slice_defining(num_vertices, subset, total_degrees):
            num_slice_defining_sets += 1
            
    # Each slice is an unordered pair {E1, E \ E1}, so we divide by 2.
    num_slices = num_slice_defining_sets // 2
    return num_slices

# --- Main analysis ---

# M(0): Based on the theorem that every cubic graph has a slice, N(G) > 0.
m_0 = "none"
print("Analysis for M(0):")
print("Every cubic graph has a slice, so N(G) is never 0. Thus, M(0) is none.")
print("-" * 20)

# M(3): Start with the smallest cubic graph, K4 (m=4).
V_k4 = 4
E_k4 = list(itertools.combinations(range(V_k4), 2))
n_k4 = count_slices(V_k4, E_k4)

print("Analysis for M(3):")
print(f"The smallest cubic graph is K4 (m=4).")
print(f"The number of slices for K4, N(K4), is {n_k4}.")
if n_k4 % 3 == 0:
    m_3 = 4
    print(f"Since {n_k4} is a multiple of 3, the smallest m is 4. M(3) = 4.")
else:
    # We would need to check m=6, but K4 works.
    m_3 = "Error in logic"
print("-" * 20)


# M(5): Check if K4 works. If not, check graphs with m=6.
print("Analysis for M(5):")
print(f"For m=4, N(K4) = {n_k4}.")
if n_k4 % 5 == 0:
    m_5 = 4
    print(f"Since {n_k4} is a multiple of 5, the smallest m is 4. M(5) = 4.")
else:
    print(f"Since {n_k4} is not a multiple of 5, we check graphs with m=6.")
    
    # Prism Graph (P6)
    V_p6 = 6
    E_p6 = [(0,1), (1,2), (2,0), (3,4), (4,5), (5,3), (0,3), (1,4), (2,5)]
    n_p6 = count_slices(V_p6, E_p6)
    print(f"For the Prism graph (m=6), N(P6) = {n_p6}.")
    
    # K_{3,3}
    V_k33 = 6
    E_k33 = []
    part1 = [0, 1, 2]
    part2 = [3, 4, 5]
    for i in part1:
        for j in part2:
            E_k33.append((i, j))
    n_k33 = count_slices(V_k33, E_k33)
    print(f"For the K_3,3 graph (m=6), N(K_3,3) = {n_k33}.")

    if n_p6 % 5 == 0 or n_k33 % 5 == 0:
        m_5 = 6
        print(f"Since N(P6)={n_p6} is a multiple of 5, a graph with m=6 exists. M(5) = 6.")
    else:
        # We would need to check m=8, but m=6 works.
        m_5 = "Error in logic"
print("-" * 20)


# Final result
print("Final Answer Summary:")
print(f"M(0) = {m_0}")
print(f"M(3) = {m_3}")
print(f"M(5) = {m_5}")
print("The final answer in the required format 'M(0),M(3),M(5)' is:")
print(f"{m_0},{m_3},{m_5}")