def solve_cubane_isomers():
    """
    This function identifies the four pairs of diagonally opposite vertices
    on the cubane structure provided in the problem image.
    """
    # Adjacency list representing the cubane graph from the product image
    adj = {
        1: [2, 4, 5],
        2: [1, 3, 6],
        3: [2, 4, 7],
        4: [1, 3, 8],
        5: [1, 6, 8],
        6: [2, 5, 7],
        7: [3, 6, 8],
        8: [4, 5, 7]
    }

    nodes = list(adj.keys())
    opposite_pairs = set()

    # Find the opposite vertex for each node (the one at distance 3)
    for start_node in nodes:
        # Use Breadth-First Search (BFS) to find the vertex at distance 3
        queue = [(start_node, 0)]  # (node, distance)
        visited = {start_node}
        
        head = 0
        while head < len(queue):
            current_node, dist = queue[head]
            head += 1

            if dist == 3:
                # Found the opposite node
                # Add the pair to the set, sorted to avoid duplicates like (1,7) and (7,1)
                pair = tuple(sorted((start_node, current_node)))
                opposite_pairs.add(pair)
                break
            
            for neighbor in adj[current_node]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append((neighbor, dist + 1))
    
    # Sort the pairs for a consistent output order
    sorted_pairs = sorted(list(opposite_pairs))

    # Format the output string as requested
    output_string = ", ".join([f"({p[0]},{p[1]})" for p in sorted_pairs])
    
    print(output_string)

solve_cubane_isomers()