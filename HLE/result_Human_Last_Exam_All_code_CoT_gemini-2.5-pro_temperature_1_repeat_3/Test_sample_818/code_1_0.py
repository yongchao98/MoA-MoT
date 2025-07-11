import itertools
from collections import deque

def get_graph_properties(adj):
    """Get number of vertices and edges from adjacency list."""
    V = len(adj)
    E = sum(len(neighbors) for neighbors in adj.values()) // 2
    return V, E

def get_incident_edges_count(adj, S):
    """Count edges incident to at least one vertex in S."""
    incident_edges = set()
    for v in S:
        for neighbor in adj[v]:
            edge = tuple(sorted((v, neighbor)))
            incident_edges.add(edge)
    return len(incident_edges)

def get_components_count(adj, S):
    """Count connected components in the subgraph induced by S."""
    if not S:
        return 0
    
    q = deque()
    visited = set()
    components = 0
    
    # Create adjacency list for the induced subgraph
    induced_adj = {v: [] for v in S}
    for u in S:
        for v in adj[u]:
            if v in S:
                induced_adj[u].append(v)

    for v in S:
        if v not in visited:
            components += 1
            q.append(v)
            visited.add(v)
            while q:
                curr = q.popleft()
                for neighbor in induced_adj[curr]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        q.append(neighbor)
    return components

def calculate_n_g(name, adj):
    """Calculate N(G) for a given graph using inclusion-exclusion."""
    V_count, E_count = get_graph_properties(adj)
    vertices = list(adj.keys())
    
    n_c_g = 0
    for i in range(V_count + 1):
        for S in itertools.combinations(vertices, i):
            S = set(S)
            s_size = len(S)
            k_s = get_components_count(adj, S)
            e_s_size = get_incident_edges_count(adj, S)
            
            term = ((-1) ** s_size) * (2 ** k_s) * (2 ** (E_count - e_s_size))
            n_c_g += term
            
    # N(G) is the number of unordered partitions, which is N_c(G) / 2
    n_g = n_c_g // 2
    # print(f"N({name}) = {n_g}")
    return n_g

def solve():
    """
    Solves for M(0), M(3), and M(5) by analyzing cubic graphs in increasing order of size.
    """
    # M(0): N(G) must be 0.
    # All cubic graphs have a 2-factor, which implies N(G) >= 1.
    # Therefore, no such graph exists.
    m0 = "none"

    # Define graphs
    K4 = {
        1: [2, 3, 4], 2: [1, 3, 4], 3: [1, 2, 4], 4: [1, 2, 3]
    }
    P6 = { # Prism graph C3 x K2
        1: [2, 3, 4], 2: [1, 3, 5], 3: [1, 2, 6],
        4: [1, 5, 6], 5: [2, 4, 6], 6: [3, 4, 5]
    }
    K3_3 = {
        1: [4, 5, 6], 2: [4, 5, 6], 3: [4, 5, 6],
        4: [1, 2, 3], 5: [1, 2, 3], 6: [1, 2, 3]
    }
    Q3 = { # Cube graph
        1: [2, 4, 5], 2: [1, 3, 6], 3: [2, 4, 7], 4: [1, 3, 8],
        5: [1, 6, 8], 6: [2, 5, 7], 7: [3, 6, 8], 8: [4, 5, 7]
    }

    # M(3)
    n_k4 = calculate_n_g("K4", K4)
    m3 = 0
    if n_k4 % 3 == 0:
        m3 = 4  # Smallest cubic graph has 4 vertices

    # M(5)
    m5 = 0
    if n_k4 % 5 == 0:
        m5 = 4
    else:
        n_p6 = calculate_n_g("P6", P6)
        n_k3_3 = calculate_n_g("K3,3", K3_3)
        if n_p6 % 5 == 0 or n_k3_3 % 5 == 0:
            m5 = 6
        else:
            n_q3 = calculate_n_g("Q3", Q3)
            if n_q3 % 5 == 0:
                m5 = 8
            else:
                m5 = "none" # Based on the problem's scope

    print(f"{m0},{m3},{m5}")

solve()