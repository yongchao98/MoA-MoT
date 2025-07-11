import collections

def solve_grid_word_search():
    """
    Finds the longest word in a grid starting with 'N' using a recursive DFS approach.
    """
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    grid_rows, grid_cols = 4, 4

    # A sample dictionary. In a real-world scenario, this would be loaded from a file.
    # It only contains words starting with 'N' since that is the problem constraint.
    word_list = ["NEON", "NOSE", "NOPES", "NEOPS", "NEOPSE"]

    # Using sets for efficient O(1) average time complexity lookups.
    word_set = set(word_list)
    prefix_set = set()
    for word in word_list:
        for i in range(1, len(word) + 1):
            prefix_set.add(word[:i])

    # This will store the results: {word: [(char, (row, col)), ...]}
    found_words = {}

    def search(row, col, path):
        """
        Recursively explores paths from a starting cell (row, col).
        'path' is a list of (character, (row, col)) tuples.
        """
        char = grid[row][col]
        current_path = path + [(char, (row, col))]
        current_word = "".join([p[0] for p in current_path])
        
        # Pruning: If the current string is not a prefix of any word, stop.
        if current_word not in prefix_set:
            return

        # If it's a valid word, record it and its path.
        if current_word in word_set:
            # We only store the first path found for a given word.
            if current_word not in found_words:
                found_words[current_word] = current_path
        
        visited_coords = {p[1] for p in current_path}

        # Explore all 8 neighbors (horizontal, vertical, diagonal).
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue  # Skip the current cell itself.
                
                new_row, new_col = row + dr, col + dc

                # Check if the neighbor is within grid bounds and not visited yet.
                if 0 <= new_row < grid_rows and 0 <= new_col < grid_cols:
                    if (new_row, new_col) not in visited_coords:
                        search(new_row, new_col, current_path)

    # Start the search from every 'N' in the grid.
    for r in range(grid_rows):
        for c in range(grid_cols):
            if grid[r][c] == 'N':
                # Start a new search with an empty initial path.
                search(r, c, [])

    if not found_words:
        print("No valid word starting with 'N' was found.")
    else:
        # Find the longest word among the ones we found.
        longest_word = max(found_words.keys(), key=len)
        path = found_words[longest_word]
        
        # Format the output as a path or "equation".
        path_str = " -> ".join([f"{char}({r},{c})" for char, (r, c) in path])
        print(f"The longest word found is '{longest_word}'.")
        print(f"Path: {path_str}")

solve_grid_word_search()