import collections

def find_para_positions_on_cubane():
    """
    Identifies the four 'para' pairs of atoms on a cubane molecule
    with a specific numbering scheme by finding vertices with a graph distance of 3.
    """
    # Define the cubane structure based on the numbering in the product image.
    # An adjacency list represents the bonds (edges) between carbon atoms (vertices).
    # Top face: 1-2-6-5, Bottom face: 4-8-7-3, Verticals: 1-4, 2-8, 5-3, 6-7
    adj = {
        1: [2, 5, 4],
        2: [1, 6, 8],
        3: [4, 7, 5],
        4: [1, 8, 3],
        5: [1, 6, 3],
        6: [2, 5, 7],
        7: [3, 8, 6],
        8: [2, 4, 7]
    }

    para_pairs = set()
    nodes = sorted(adj.keys())

    # Find all pairs of nodes with a shortest path of 3
    for start_node in nodes:
        # BFS to find distances from start_node
        q = collections.deque([(start_node, 0)]) # (node, distance)
        visited = {start_node}
        distances = {start_node: 0}

        while q:
            current_node, distance = q.popleft()

            # The body diagonal of a cube has length 3
            if distance == 3:
                # We found a para pair. Add it to the set.
                # Sorting ensures (a,b) is the same as (b,a)
                pair = tuple(sorted((start_node, current_node)))
                para_pairs.add(pair)
                # No need to search further from this path
                continue

            for neighbor in adj[current_node]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    distances[neighbor] = distance + 1
                    q.append((neighbor, distance + 1))
    
    # Sort the pairs for a consistent output order
    sorted_pairs = sorted(list(para_pairs))

    # Format the output as requested
    output_string = ", ".join([f"({p[0]},{p[1]})" for p in sorted_pairs])
    print("The four theoretically possible pairs of substituted carbon atoms are:")
    print(output_string)


find_para_positions_on_cubane()
<<< (1,7), (2,3), (4,6), (5,8) >>>