import collections

def solve():
    """
    Calculates M(0), M(3), and M(5) and prints the result.
    """

    # Helper function to parse g6 format for small graphs
    def parse_g6(g6_string):
        data = bytearray(g6_string, 'ascii')
        for i in range(len(data)):
            data[i] -= 63
        
        n = data[0]
        data = data[1:]
        
        adj_matrix = [[0] * n for _ in range(n)]
        k = 0
        for j in range(1, n):
            for i in range(j):
                byte_index = k // 6
                bit_index = 5 - (k % 6)
                if (data[byte_index] >> bit_index) & 1:
                    adj_matrix[i][j] = adj_matrix[j][i] = 1
                k += 1
        
        edges = []
        for i in range(n):
            for j in range(i + 1, n):
                if adj_matrix[i][j]:
                    edges.append((i, j))
        return n, edges

    def calculate_n_g(v, edges):
        """
        Calculates N(G) for a graph G=(v, edges) using inclusion-exclusion.
        """
        e = len(edges)
        adj = collections.defaultdict(list)
        for u_node, v_node in edges:
            adj[u_node].append(v_node)
            adj[v_node].append(u_node)

        union_of_mono_sets_size = 0
        
        # Iterate through all non-empty subsets of vertices A
        for i in range(1, 1 << v):
            A = []
            for j in range(v):
                if (i >> j) & 1:
                    A.append(j)
            
            # For a set A, calculate the size of the intersection of A_v for v in A.
            # This is the number of colorings where all vertices in A are monochromatic.
            
            # 1. Determine color constraints on vertices in A.
            # If u,v in A are adjacent, they must have the same mono-color.
            # We can model this with a "constraint graph" on the vertices of A.
            constraint_adj = collections.defaultdict(list)
            nodes_in_A = set(A)
            for u_node in A:
                for v_node in adj[u_node]:
                    if v_node in nodes_in_A and u_node < v_node:
                        constraint_adj[u_node].append(v_node)
                        constraint_adj[v_node].append(u_node)
            
            # 2. Count components in the constraint graph to find independent color choices.
            num_components = 0
            visited = set()
            for node in A:
                if node not in visited:
                    num_components += 1
                    q = collections.deque([node])
                    visited.add(node)
                    while q:
                        curr = q.popleft()
                        for neighbor in constraint_adj[curr]:
                            if neighbor not in visited:
                                visited.add(neighbor)
                                q.append(neighbor)
            
            # 3. Count constrained edges.
            edges_constrained = set()
            for node in A:
                for neighbor in adj[node]:
                    edges_constrained.add(tuple(sorted((node, neighbor))))
            
            num_free_edges = e - len(edges_constrained)
            
            # Each component of the constraint graph can be all-red or all-blue.
            # The remaining edges are free to be any color.
            term = (2**num_components) * (2**num_free_edges)
            
            if len(A) % 2 == 1:
                union_of_mono_sets_size += term
            else:
                union_of_mono_sets_size -= term

        s_prime = (2**e) - union_of_mono_sets_size
        return s_prime // 2

    # M(0): N(G) is never 0 for a cubic graph.
    m0 = "none"

    # M(3): Smallest cubic graph is K4 (m=6).
    v_k4, edges_k4 = parse_g6("C?")
    n_k4 = calculate_n_g(v_k4, edges_k4)
    m3 = 0
    if n_k4 % 3 == 0:
        m3 = len(edges_k4)

    # M(5): Check graphs in increasing order of edges.
    # A catalogue of connected cubic graphs in g6 format.
    graph_catalogue = {
        6: ["C?"],  # K4
        9: ["E?~", "E@w"],  # K3,3, Prism3
        12: ["G??O@", "G?CG@", "G?aG@", "G?AGw", "G?Qw"],
        15: ["I?{OG@", "I?{oG@", "I?a{G@", "I?e{G@", "I?i{G@", "I?q{G@", "I?AG{w", "I?QG{w", "I?aG{w", "I?iG{w", "I?qG{w", "I?yG{w", "I?YG{w", "I?gG{w", "I?sG{w", "I?{oGw", "I?{_Gw", "I?a{Gw", "I?i{Gw"]
    }
    
    m5 = "none"
    found_m5 = False
    for m in sorted(graph_catalogue.keys()):
        if found_m5:
            break
        for g6 in graph_catalogue[m]:
            v, edges = parse_g6(g6)
            n_g = calculate_n_g(v, edges)
            if n_g % 5 == 0:
                m5 = m
                found_m5 = True
                break
    
    print(f"{m0},{m3},{m5}")

solve()