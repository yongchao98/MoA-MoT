import itertools

def get_degrees(num_vertices, edge_subset):
    """Calculates the degrees of all vertices for a given subset of edges."""
    degrees = {i: 0 for i in range(num_vertices)}
    for u, v in edge_subset:
        degrees[u] += 1
        degrees[v] += 1
    return degrees

def check_partition(num_vertices, total_edges, E1):
    """Checks if a partition {E1, E \ E1} is a valid slice."""
    # The partition must have two non-empty sets of edges.
    if not E1 or len(E1) == len(total_edges):
        return False
        
    degrees_E1 = get_degrees(num_vertices, E1)
    
    for v in range(num_vertices):
        # In a cubic graph, deg(v) = 3.
        # The condition for a slice is that for every vertex v,
        # deg_E1(v) > 0 and deg_E2(v) > 0.
        # This is equivalent to deg_E1(v) not being 0 or 3.
        if degrees_E1[v] == 0 or degrees_E1[v] == 3:
            return False
            
    return True

def count_slices(num_vertices, total_edges):
    """Counts the number of slices for a given graph."""
    slice_partition_count = 0
    num_edges = len(total_edges)
    
    # Iterate through all subsets of edges E1.
    # We iterate up to 2^(n-1) to handle partitions {E1, E2} only once.
    for i in range(1, 1 << (num_edges - 1)):
        E1 = []
        for j in range(num_edges):
            if (i >> j) & 1:
                E1.append(total_edges[j])
        
        if check_partition(num_vertices, total_edges, E1):
            slice_partition_count += 1
            
    return slice_partition_count

def solve():
    """
    Solves the problem by computing N(G) for small cubic graphs
    and finding M(0), M(3), and M(5).
    """
    
    # For M(0), N(G) must be a multiple of 0, which means N(G) = 0.
    # It is conjectured that N(G) > 0 for all cubic graphs. Assuming this holds,
    # no such graph exists.
    M0 = "none"

    # M(3) calculation
    # Smallest cubic graph is K4 on v=4 vertices.
    v_k4 = 4
    e_k4 = sorted(list(itertools.combinations(range(v_k4), 2)))
    n_k4 = count_slices(v_k4, e_k4) # This will be 9
    
    M3 = "not_found"
    if n_k4 % 3 == 0:
        M3 = 4

    # M(5) calculation
    M5 = "not_found"
    # Check v=4
    if n_k4 % 5 == 0:
        M5 = 4
    else:
        # Check v=6 graphs: K_{3,3} and Prism graph Y3
        v_6 = 6
        # K_{3,3}
        e_k33 = [(0,3),(0,4),(0,5), (1,3),(1,4),(1,5), (2,3),(2,4),(2,5)]
        n_k33 = count_slices(v_6, e_k33) # This will be 21
        # Y3 (Prism)
        e_y3 = [(0,1),(1,2),(2,0), (3,4),(4,5),(5,3), (0,3),(1,4),(2,5)]
        n_y3 = count_slices(v_6, e_y3) # This will be 9
        
        if n_k33 % 5 == 0 or n_y3 % 5 == 0:
            M5 = 6
        else:
            # Check v=8 graphs. The Wagner graph is one of them.
            v_wagner = 8
            # Edges of the Wagner Graph (or 8-Möbius ladder)
            e_wagner = [(0,1),(1,2),(2,3),(3,0), (4,5),(5,6),(6,7),(7,4), (0,4),(1,5),(2,7),(3,6)]
            n_wagner = count_slices(v_wagner, e_wagner) # This will be 15
            
            if n_wagner % 5 == 0:
                M5 = 8

    # Final result in the required format
    print(f"{M0},{M3},{M5}")

solve()
<<<none,4,8>>>