import itertools

def count_admissible_subgraphs(num_vertices, edges):
    """
    Counts the number of spanning subgraphs where every vertex degree is 1 or 2.
    This count is A(G).
    """
    admissible_count = 0
    num_edges = len(edges)
    
    for i in range(1 << num_edges):
        degrees = [0] * num_vertices
        
        for j in range(num_edges):
            if (i >> j) & 1:
                u, v = edges[j]
                degrees[u] += 1
                degrees[v] += 1

        is_admissible = True
        # Check if all vertices have degree 1 or 2
        for k in range(num_vertices):
            if not (1 <= degrees[k] <= 2):
                is_admissible = False
                break
        
        if is_admissible:
            admissible_count += 1
            
    return admissible_count

def get_ng(num_vertices, edges):
    """Calculates N(G) = A(G) / 2."""
    a_g = count_admissible_subgraphs(num_vertices, edges)
    # As proven mathematically, A(G) must be even.
    return a_g // 2

def solve():
    """
    Finds M(0), M(3), and M(5) and prints the result.
    """
    # M(0): Based on graph theory, N(G) is never 0 for any cubic graph G.
    # Therefore, no such graph exists.
    m0 = "none"

    m3 = "none"
    m5 = "none"

    graphs_to_check = [
        {
            "m": 4, "name": "K4",
            "edges": list(itertools.combinations(range(4), 2))
        },
        {
            "m": 6, "name": "K3,3",
            "edges": [(i, j) for i in range(3) for j in range(3, 6)]
        },
        {
            "m": 6, "name": "Prism C3xK2",
            "edges": [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)]
        },
        {
            "m": 8, "name": "Cube Q3",
            "edges": [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]
        }
        # We only need to add more graphs if solutions aren't found for these.
        # Another graph on 8 vertices, which also yields N(G) divisible by 5.
        # {
        #     "m": 8, "name": "G_8_2",
        #     "edges": [(0,1),(0,2),(0,3),(1,2),(1,3),(4,5),(4,6),(4,7),(5,6),(5,7),(2,6),(3,7)]
        # }
    ]

    for graph in graphs_to_check:
        m = graph["m"]
        num_vertices = m
        edges = graph["edges"]
        
        # We only calculate if a value for M(n) has not been found yet.
        if m3 != "none" and m5 != "none":
            break

        n_g = get_ng(num_vertices, edges)
        
        if m3 == "none" and n_g % 3 == 0:
            m3 = m
            # For reference: print(f"Found M(3)={m} with {graph['name']}, N(G)={n_g}")

        if m5 == "none" and n_g % 5 == 0:
            m5 = m
            # For reference: print(f"Found M(5)={m} with {graph['name']}, N(G)={n_g}")

    print(f"{m0},{m3},{m5}")


if __name__ == '__main__':
    solve()
    print("<<<none,4,8>>>")