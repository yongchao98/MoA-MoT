import collections

def count_arrangements():
    """
    Calculates the number of valid arrangements for the letters L, N, S, and W
    by finding all Hamiltonian paths in a directed graph derived from connection rules.
    """
    # The graph of valid connections is defined as a dictionary.
    # This graph is derived from an intuitive model of how letters connect,
    # patched with the specific rules given in the problem's example ("S" -> "W" is valid, "W" -> "S" is not).
    # L(ends B)->N(starts B) -> OK
    # S(ends B)->N(starts B) -> OK
    # N(ends T)->L(starts T),S(starts T),W(starts T) -> OK
    # W(ends T)->L(starts T) -> OK
    # The example S->W is added and W->S is removed from the intuitive model.
    graph = {
        'L': ['N'],
        'N': ['L', 'S', 'W'],
        'S': ['N', 'W'],  # 'W' is added based on the example.
        'W': ['L']        # 'S' is removed based on the example.
    }

    letters = ['L', 'N', 'S', 'W']
    all_valid_paths = []

    def find_paths_recursively(current_path):
        """
        Uses recursion to explore the graph and find all paths that use each letter exactly once.
        """
        # If the path includes all 4 letters, it's a valid arrangement.
        if len(current_path) == 4:
            all_valid_paths.append(current_path)
            return

        last_letter_in_path = current_path[-1]

        # Explore next possible letters from the last one in the path.
        if last_letter_in_path in graph:
            for neighbor in graph[last_letter_in_path]:
                if neighbor not in current_path:  # Ensure we don't visit the same letter twice.
                    find_paths_recursively(current_path + [neighbor])

    # Start the search from each of the four possible starting letters.
    for start_letter in letters:
        find_paths_recursively([start_letter])
    
    # Deconstruct the total number of paths into a sum, as requested.
    paths_by_start_node = collections.defaultdict(int)
    for path in all_valid_paths:
        paths_by_start_node[path[0]] += 1

    # Print the "equation" showing the number of paths starting with each letter.
    # The order is fixed for consistent output.
    l_paths = paths_by_start_node.get('L', 0)
    n_paths = paths_by_start_node.get('N', 0)
    s_paths = paths_by_start_node.get('S', 0)
    w_paths = paths_by_start_node.get('W', 0)
    
    total_paths = len(all_valid_paths)

    # Example valid paths found: [W, L, N, S], [S, N, W, L], etc.
    # Printing the count as a sum for clarity.
    print(f"{l_paths} + {n_paths} + {s_paths} + {w_paths} = {total_paths}")


count_arrangements()
<<<5>>>