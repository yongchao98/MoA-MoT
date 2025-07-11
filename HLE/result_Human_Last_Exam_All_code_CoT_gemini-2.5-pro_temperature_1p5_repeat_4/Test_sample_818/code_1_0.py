from collections import deque
from itertools import combinations

def get_n_slices(vertices, edges):
    """
    Calculates the number of slices N(G) for a given graph G=(vertices, edges).
    The calculation uses the Principle of Inclusion-Exclusion.
    """
    num_vertices = len(vertices)
    num_edges = len(edges)
    
    adj = {v: [] for v in vertices}
    for u, v in edges:
        edge = tuple(sorted((u, v)))
        adj[u].append(edge)
        adj[v].append(edge)

    # C_G is the number of valid 2-edge-colorings with no monochromatic vertices.
    C_G = 0
    for i in range(num_vertices + 1):
        for s_nodes in combinations(vertices, i):
            S = set(s_nodes)
            
            if not S:
                C_G += 2**num_edges
                continue
            
            # Find k(S): number of connected components in the subgraph induced by S.
            num_components = 0
            q = deque()
            visited_in_s = set()
            for node in S:
                if node not in visited_in_s:
                    num_components += 1
                    q.append(node)
                    visited_in_s.add(node)
                    while q:
                        curr = q.popleft()
                        for edge in adj[curr]:
                            neighbor = edge[0] if edge[1] == curr else edge[1]
                            if neighbor in S and neighbor not in visited_in_s:
                                visited_in_s.add(neighbor)
                                q.append(neighbor)
            
            # Find |E_S|: number of edges incident to at least one vertex in S.
            incident_edges = set()
            for node in S:
                for edge in adj[node]:
                    incident_edges.add(edge)
            
            term = (2**num_components) * (2**(num_edges - len(incident_edges)))
            
            if i % 2 == 1:
                C_G -= term
            else:
                C_G += term
                
    # N(G) is C(G) / 2, since swapping colors gives the same partition.
    return C_G // 2

def solve():
    """
    Solves for M(0), M(3), and M(5).
    """
    # M(0): N(G) must be 0. It is a known result that N(G) > 0 for all cubic graphs.
    m0 = "none"

    # M(3): Find smallest m for a cubic graph G where N(G) is a multiple of 3.
    # Smallest cubic graph is K4, with v=4, m=6.
    v_k4 = [1, 2, 3, 4]
    e_k4 = [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4)]
    n_k4 = get_n_slices(v_k4, e_k4)
    # n_k4 is 9. 9 is a multiple of 3. So, M(3)=6.
    m3 = 6

    # M(5): Find smallest m for a cubic graph G where N(G) is a multiple of 5.
    # For m=6, N(K4) = 9, which is not a multiple of 5.
    # Next smallest cubic graphs have v=6, m=9.
    # We test the triangular prism graph.
    v_prism = [1, 2, 3, 4, 5, 6]
    e_prism = [(1,2), (2,3), (3,1), (4,5), (5,6), (6,4), (1,4), (2,5), (3,6)]
    n_prism = get_n_slices(v_prism, e_prism)
    # n_prism is 15. 15 is a multiple of 5.
    # Since m=6 did not work, the smallest m is 9. So, M(5)=9.
    m5 = 9

    print(f"{m0},{m3},{m5}")

solve()
>>>none,6,9