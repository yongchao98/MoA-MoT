import collections

def find_para_pairs():
    """
    Identifies the para-substituted pairs on the cubane product.

    The cubane's connectivity is defined as per the product image.
    In a cube, para positions are diagonally opposite vertices (body diagonals).
    The shortest path between them along the edges is 3.
    """

    # Define the connections (edges) of the cube based on the product image
    edges = {
        1: [2, 4, 5],
        2: [1, 3, 6],
        3: [2, 4, 7],
        4: [1, 3, 8],
        5: [1, 6, 8],
        6: [2, 5, 7],
        7: [3, 6, 8],
        8: [4, 5, 7]
    }

    vertices = list(edges.keys())
    para_pairs = []
    processed_vertices = set()

    for start_node in vertices:
        if start_node in processed_vertices:
            continue

        # Use Breadth-First Search (BFS) to find distances from start_node
        queue = collections.deque([(start_node, [start_node])])
        visited = {start_node}
        
        para_partner = None
        while queue:
            current_node, path = queue.popleft()
            
            # If path length is 4 (distance is 3), we found the para partner
            if len(path) == 4:
                para_partner = current_node
                break

            for neighbor in edges[current_node]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    new_path = list(path)
                    new_path.append(neighbor)
                    queue.append((neighbor, new_path))
        
        if para_partner:
            pair = tuple(sorted((start_node, para_partner)))
            para_pairs.append(pair)
            processed_vertices.add(start_node)
            processed_vertices.add(para_partner)

    # Format the output string as requested
    output_string = ", ".join([f"({p[0]},{p[1]})" for p in sorted(para_pairs)])
    print(output_string)

# Run the function to get the answer
find_para_pairs()