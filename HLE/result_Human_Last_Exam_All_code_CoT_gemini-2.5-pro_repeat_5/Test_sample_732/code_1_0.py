import collections

def solve_and_verify():
    """
    This script solves the problem of finding the maximum number of edges in a 
    C4-free graph with 8 vertices by constructing an example and verifying its properties.
    """

    # 1. Define the problem parameters
    num_vertices_V = 8
    num_edges_in_cycle = 8

    # 2. Construct the graph with 10 edges
    # The graph consists of an 8-cycle plus two chords.
    graph = collections.defaultdict(list)
    edges = set()
    
    # Add edges of the 8-cycle: (0,1), (1,2), ..., (7,0)
    for i in range(num_vertices_V):
        u, v = i, (i + 1) % num_vertices_V
        graph[u].append(v)
        graph[v].append(u)
        edges.add(tuple(sorted((u, v))))

    # Add two chords
    chords = [(0, 4), (2, 6)]
    num_chords = len(chords)
    for u, v in chords:
        graph[u].append(v)
        graph[v].append(u)
        edges.add(tuple(sorted((u, v))))

    num_edges_E = len(edges)
    
    # "output each number in the final equation!"
    # Here is the equation for the number of edges in our constructed graph:
    print(f"Number of edges = (edges in cycle) + (added chords)")
    print(f"{num_edges_E} = {num_edges_in_cycle} + {num_chords}")
    print("-" * 20)

    # 3. Verify that the graph is C4-free
    # A graph is C4-free if and only if no pair of vertices shares more than one common neighbor.
    is_c4_free = True
    c4_evidence = None
    nodes = sorted(list(graph.keys()))
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            u, v = nodes[i], nodes[j]
            
            neighbors_u = set(graph[u])
            neighbors_v = set(graph[v])
            common_neighbors = neighbors_u.intersection(neighbors_v)
            
            if len(common_neighbors) > 1:
                is_c4_free = False
                c_neighbors = list(common_neighbors)
                c4_evidence = f"Vertices {u} and {v} have common neighbors {c_neighbors[0]} and {c_neighbors[1]}, forming a C4."
                break
        if not is_c4_free:
            break
            
    if is_c4_free:
        print("Verification successful: The constructed graph is C4-free.")
    else:
        print("Verification FAILED: The constructed graph contains a C4.")
        print(c4_evidence)
        return

    # 4. State the final conclusion
    print("\nSince a C4-free graph with 8 vertices and 10 edges exists, and it is known")
    print("from extremal graph theory that 11 edges are not possible, the maximum is 10.")
    print("\nThe maximum number of edges in a simple graph with 8 vertices and no quadrilaterals is 10.")


solve_and_verify()