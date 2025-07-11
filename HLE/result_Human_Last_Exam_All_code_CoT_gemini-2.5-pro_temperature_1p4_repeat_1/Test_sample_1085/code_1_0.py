import collections

def find_prime_paths():
    """
    This function searches for all "Prime Paths" of exactly 4 moves
    from (1, 1) to (5, 7) in the PrimeGrid+1 space.
    """
    
    # The ordered set of valid indices for coordinates, based on the problem's scope.
    # It includes 1 and the prime numbers relevant to the start and end points.
    indices = [1, 2, 3, 5, 7]
    index_map = {val: i for i, val in enumerate(indices)}

    def get_neighbors(node):
        """Calculates the adjacent Prime Intersections for a given node."""
        x, y = node
        neighbors = []
        
        # This check is for safety, though all generated nodes should be valid.
        if x not in index_map or y not in index_map:
            return []

        x_idx = index_map[x]
        y_idx = index_map[y]

        # Horizontal neighbors (previous and next in the prime sequence)
        if x_idx > 0:
            neighbors.append((indices[x_idx - 1], y))
        if x_idx < len(indices) - 1:
            neighbors.append((indices[x_idx + 1], y))

        # Vertical neighbors (previous and next in the prime sequence)
        if y_idx > 0:
            neighbors.append((x, indices[y_idx - 1]))
        if y_idx < len(indices) - 1:
            neighbors.append((x, indices[y_idx + 1]))
            
        return neighbors

    start_node = (1, 1)
    end_node = (5, 7)
    num_moves = 4
    
    # We use a breadth-first search (BFS) to explore all paths of length 4.
    # The queue will store paths, where each path is a list of nodes.
    queue = collections.deque([[start_node]])
    
    found_paths = []
    
    # The search proceeds level by level, where a level corresponds to one move.
    while queue:
        path = queue.popleft()
        
        # A path with 'num_moves' has 'num_moves + 1' nodes.
        # If the path has reached the desired length, we check its destination.
        if len(path) == num_moves + 1:
            if path[-1] == end_node:
                found_paths.append(path)
            # Stop exploring this branch as we must have exactly 4 moves.
            continue

        # If the path is still shorter than 4 moves, extend it.
        current_node = path[-1]
        for neighbor in get_neighbors(current_node):
            new_path = list(path)
            new_path.append(neighbor)
            queue.append(new_path)
            
    # Output the results.
    print(f"Searching for paths from {start_node} to {end_node} with exactly {num_moves} moves.\n")
    
    if not found_paths:
        print("No paths were found that match the criteria.")
        # The "equation" for zero paths is trivial.
        print("\nFinal Equation:")
        print("0 = 0")
    else:
        print(f"Found {len(found_paths)} distinct paths:")
        for i, path in enumerate(found_paths):
            path_str = " -> ".join(map(str, path))
            print(f"Path {i+1}: {path_str}")
        
        # Construct the final equation string: 1 + 1 + ... = Total
        equation_str = " + ".join(["1"] * len(found_paths))
        print("\nFinal Equation:")
        print(f"{equation_str} = {len(found_paths)}")

    print(f"\nTotal number of distinct Prime Paths: {len(found_paths)}")


# Execute the function to find and print the number of paths.
find_prime_paths()
<<<0>>>