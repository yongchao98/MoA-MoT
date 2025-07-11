import collections

def solve_prime_path():
    """
    Calculates the number of distinct Prime Paths from (1,1) to (5,7) in 4 moves.
    """
    # Step 1 & 3: Define the grid indices and the adjacency rules.
    # Adjacency based on the interpretation that 1 is adjacent to all primes,
    # and primes are adjacent to their prime neighbors and 1.
    adj = {
        1: {2, 3, 5, 7},
        2: {1, 3},
        3: {1, 2, 5},
        5: {1, 3, 7},
        7: {1, 5}
    }

    start_node = (1, 1)
    end_node = (5, 7)
    path_length = 4
    
    all_paths = []
    
    # Step 4: Use a recursive DFS to find all paths of the specified length.
    # A queue for the DFS stores tuples of (current_path_list).
    queue = collections.deque([[start_node]])

    while queue:
        path = queue.popleft()
        
        if len(path) == path_length + 1:
            if path[-1] == end_node:
                all_paths.append(tuple(path))
            continue

        last_node = path[-1]
        x, y = last_node

        # Find horizontal neighbors
        for nx in adj.get(x, []):
            neighbor = (nx, y)
            new_path = list(path)
            new_path.append(neighbor)
            queue.append(new_path)
            
        # Find vertical neighbors
        for ny in adj.get(y, []):
            neighbor = (x, ny)
            new_path = list(path)
            new_path.append(neighbor)
            queue.append(new_path)

    # Step 5: Categorize the found paths.
    # Using sets to ensure we count distinct path sequences.
    unique_paths = set(all_paths)
    
    cat_c = 0  # No backtracking: P0, P1, P2, P3, P4 are all unique
    cat_a_only = 0 # Backtracking P0 -> P1 -> P0 -> ...
    cat_b_only = 0 # Backtracking P0 -> P1 -> P2 -> P1 -> ...
    cat_ab_overlap = 0 # Overlap: P0 -> P1 -> P0 -> P1 -> ...

    for p in unique_paths:
        is_a = (p[0] == p[2])
        is_b = (p[1] == p[3])
        
        # Path has no repeated nodes
        if len(set(p)) == 5:
            cat_c += 1
        elif is_a and is_b:
            cat_ab_overlap += 1
        elif is_a:
            cat_a_only += 1
        elif is_b:
            cat_b_only += 1
            
    total_paths = len(unique_paths)

    # Step 6: Print the results, showing the "final equation".
    print("Analyzing distinct paths of 4 moves from (1, 1) to (5, 7):")
    print("-" * 60)
    print(f"Paths with no backtracking (e.g., (1,1)->(1,3)->(1,5)->(1,7)->(5,7)): {cat_c}")
    print(f"Paths with P0->P1->P0... backtracking (A only): {cat_a_only}")
    print(f"Paths with ...P1->P2->P1... backtracking (B only): {cat_b_only}")
    print(f"Paths with P0->P1->P0->P1... backtracking (A and B overlap): {cat_ab_overlap}")
    print("-" * 60)
    print("The final equation for the total number of distinct paths is:")
    print(f"{cat_c} + {cat_a_only} + {cat_b_only} + {cat_ab_overlap} = {total_paths}")
    print(f"Total number of distinct Prime Paths: {total_paths}")

solve_prime_path()