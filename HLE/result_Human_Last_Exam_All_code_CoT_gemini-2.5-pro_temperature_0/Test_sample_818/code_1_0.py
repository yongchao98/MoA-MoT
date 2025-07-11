import itertools

def count_slices(num_vertices, edges):
    """
    Counts the number of slices in a graph by checking all edge partitions.
    A slice is a partition of edges into two classes (C1, C2) such that
    every vertex is incident to at least one edge in each class.
    """
    num_edges = len(edges)
    
    # We iterate through all 2^num_edges possible colorings of the edges.
    # A coloring is an assignment of each edge to one of two classes.
    valid_colorings = 0
    for i in range(2**num_edges):
        # Generate the i-th coloring
        coloring = {}
        temp_i = i
        for j in range(num_edges):
            coloring[edges[j]] = temp_i % 2
            temp_i //= 2
            
        is_valid_coloring = True
        for v in range(num_vertices):
            incident_colors = set()
            for edge in edges:
                if v in edge:
                    incident_colors.add(coloring[edge])
            
            # For a slice, each vertex must touch edges of both colors.
            if len(incident_colors) < 2:
                is_valid_coloring = False
                break
        
        if is_valid_coloring:
            valid_colorings += 1
            
    # Each slice corresponds to two complementary colorings (e.g., all 0s become 1s and vice-versa).
    # So we divide the number of valid colorings by 2.
    return valid_colorings // 2

def solve():
    """
    Determines M(0), M(3), and M(5) by analyzing small cubic graphs.
    """
    # M(0): A multiple of 0 must be 0. It can be shown that any cubic graph G
    # has N(G) >= 1. Thus, no graph G has N(G) being a multiple of 0.
    m0 = "none"

    # M(3): Find smallest m where N(G) is a multiple of 3.
    # Start with the smallest cubic graph, K4 (m=4).
    v_k4 = 4
    e_k4 = [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]
    n_k4 = count_slices(v_k4, e_k4) # Expected: 9
    
    m3 = "none"
    if n_k4 % 3 == 0:
        m3 = 4

    # M(5): Find smallest m where N(G) is a multiple of 5.
    # Check m=4 first.
    m5 = "none"
    if n_k4 % 5 == 0:
        m5 = 4
    else:
        # Check m=6 graphs. There are two: Prism graph and K_{3,3}.
        # We only need to find one that satisfies the condition.
        
        # K_{3,3} graph
        v_k33 = 6
        e_k33 = [(0,3), (0,4), (0,5), (1,3), (1,4), (1,5), (2,3), (2,4), (2,5)]
        n_k33 = count_slices(v_k33, e_k33) # Expected: 30
        
        # Prism graph (C3 x K2)
        v_p6 = 6
        e_p6 = [(0,1),(1,2),(2,0), (3,4),(4,5),(5,3), (0,3),(1,4),(2,5)]
        n_p6 = count_slices(v_p6, e_p6) # Expected: 21

        if n_k33 % 5 == 0 or n_p6 % 5 == 0:
            m5 = 6
            
    print(f"{m0},{m3},{m5}")

solve()
<<<none,4,6>>>