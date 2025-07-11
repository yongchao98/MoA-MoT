def find_prime_paths():
    """
    Calculates the number of distinct 4-move paths from (1, 1) to (5, 7)
    in the PrimeGrid+1 space.
    """
    
    # The set of valid coordinates includes 1 and prime numbers.
    # For the scope of this problem, we only need primes up to 7.
    prime_coords = [1, 2, 3, 5, 7]
    
    start_node = (1, 1)
    end_node = (5, 7)
    required_moves = 4
    
    # This list will store any valid paths that are found.
    # A path is a sequence of intersections (nodes).
    found_paths = []

    def get_neighbors(node):
        """
        Finds all adjacent Prime Intersections to a given node.
        Adjacency is defined by the next/previous coordinate in the prime_coords list.
        """
        x, y = node
        neighbors = []
        
        try:
            x_idx = prime_coords.index(x)
            y_idx = prime_coords.index(y)
        except ValueError:
            # This should not be reached with a valid path.
            return []

        # Find vertical neighbors (change in x)
        if x_idx > 0:
            neighbors.append((prime_coords[x_idx - 1], y))
        if x_idx < len(prime_coords) - 1:
            neighbors.append((prime_coords[x_idx + 1], y))
            
        # Find horizontal neighbors (change in y)
        if y_idx > 0:
            neighbors.append((x, prime_coords[y_idx - 1]))
        if y_idx < len(prime_coords) - 1:
            neighbors.append((x, prime_coords[y_idx + 1]))
            
        return neighbors

    def search_paths_recursive(current_path):
        """
        A recursive function to explore all paths of a specific length.
        """
        # If the path has the required number of moves (path length = moves + 1),
        # check if it ends at the target destination.
        if len(current_path) - 1 == required_moves:
            if current_path[-1] == end_node:
                found_paths.append(current_path)
            return

        # If the path is already too long, stop exploring.
        if len(current_path) - 1 > required_moves:
            return

        # Explore from the last node in the current path.
        last_node = current_path[-1]
        for neighbor in get_neighbors(last_node):
            new_path = current_path + [neighbor]
            search_paths_recursive(new_path)

    # Start the search from the beginning node.
    search_paths_recursive([start_node])
    
    # Output the results as requested.
    count = len(found_paths)
    
    if count == 0:
        print("No paths of length 4 were found from (1, 1) to (5, 7).")
        print("Final Equation: 0 = 0")
    else:
        # This part of the code would run if paths were found.
        # It generates an equation like "1 + 1 + ... = count".
        for i, path in enumerate(found_paths):
            path_str = " -> ".join(map(str, path))
            print(f"Path {i+1}: {path_str}")
        
        equation_parts = ["1"] * count
        print(f"Final Equation: {' + '.join(equation_parts)} = {count}")
    
    print(f"\nThe total number of distinct Prime Paths is: {count}")

# Execute the function
find_prime_paths()
<<<0>>>