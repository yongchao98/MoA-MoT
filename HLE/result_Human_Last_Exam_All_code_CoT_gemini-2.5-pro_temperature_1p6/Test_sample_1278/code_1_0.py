import collections

def count_cycles():
    """
    This function defines a specific planar graph with 9 vertices and 16 edges
    and counts all simple cycles of length >= 3 within it.

    The graph is a 3x3 grid with four additional diagonal edges chosen to
    maximize the number of triangular faces, which is a heuristic for maximizing
    the total number of cycles.
    """
    # Adjacency list for the graph. Vertices are numbered 0 to 8,
    # corresponding to a 3x3 grid layout:
    # 0 1 2
    # 3 4 5
    # 6 7 8
    adj = collections.defaultdict(list)
    def add_edge(u, v):
        adj[u].append(v)
        adj[v].append(u)

    # Base 3x3 grid edges (12 edges)
    add_edge(0, 1); add_edge(1, 2)
    add_edge(3, 4); add_edge(4, 5)
    add_edge(6, 7); add_edge(7, 8)
    add_edge(0, 3); add_edge(3, 6)
    add_edge(1, 4); add_edge(4, 7)
    add_edge(2, 5); add_edge(5, 8)

    # Diagonal edges added to the 4 inner squares (4 edges)
    add_edge(0, 4)  # Corresponds to v11-v22
    add_edge(1, 5)  # Corresponds to v12-v23
    add_edge(3, 7)  # Corresponds to v21-v32
    add_edge(4, 8)  # Corresponds to v22-v33

    # Sort adjacency lists to have a consistent search order
    for i in range(9):
        adj[i].sort()

    cycles_by_length = collections.defaultdict(int)
    
    # Backtracking function to find cycles
    # To avoid duplicates, we enforce that cycles are explored from their
    # lowest-indexed node, and the path only includes higher-indexed nodes.
    def find_cycles_from_node(path):
        start_node = path[0]
        current_node = path[-1]
        
        for neighbor in adj[current_node]:
            # If the neighbor closes the cycle and it's long enough
            if neighbor == start_node and len(path) >= 3:
                cycle_len = len(path)
                cycles_by_length[cycle_len] += 1
            # Continue the path with nodes greater than the start node
            elif neighbor > start_node and neighbor not in path:
                find_cycles_from_node(path + [neighbor])

    # Start the search from each node
    for i in range(9):
        find_cycles_from_node([i])

    total_cycles = sum(cycles_by_length.values())
    
    equation_parts = []
    for length in sorted(cycles_by_length.keys()):
        count = cycles_by_length[length]
        equation_parts.append(f"{count}")

    print("Equation for the number of standoffs (cycles):")
    # Output each number in the final equation as requested
    equation_string = " + ".join(equation_parts)
    print(f"Counts per cycle length (C3 to C8): {equation_string}")
    print(f"Total number of Mexican standoffs: {total_cycles}")


count_cycles()