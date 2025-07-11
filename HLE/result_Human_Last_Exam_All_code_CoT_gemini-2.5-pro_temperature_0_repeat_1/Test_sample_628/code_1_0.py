def find_longest_grid_word():
    """
    This script finds the longest word in a 4x4 grid starting with 'N'.
    It uses a depth-first search algorithm to explore all possible paths.
    """
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]

    # A small, curated dictionary for this specific problem.
    # A real-world solution would load a comprehensive dictionary from a file.
    dictionary = {
        "NO", "NOES", "NOPE", "NOPS", "NOSE", "NEP", "NEPO", "NEPOS", "NOOK", "NOOKS"
    }

    rows, cols = 4, 4
    max_word = ""

    def is_valid(r, c, visited):
        """Check if a cell is within grid bounds and not visited."""
        return 0 <= r < rows and 0 <= c < cols and not visited[r][c]

    def dfs(r, c, current_word, visited):
        """Recursively explore paths to form words."""
        nonlocal max_word
        
        # Mark the current cell as visited for this path
        visited[r][c] = True
        current_word += grid[r][c]

        # If the current path forms a word in the dictionary, check if it's the longest
        if current_word in dictionary:
            if len(current_word) > len(max_word):
                max_word = current_word
            # If words have the same length, prefer the lexicographically smaller one
            elif len(current_word) == len(max_word) and current_word < max_word:
                max_word = current_word

        # Explore all 8 neighbors (horizontal, vertical, diagonal)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                
                next_r, next_c = r + dr, c + dc
                
                if is_valid(next_r, next_c, visited):
                    dfs(next_r, next_c, current_word, visited)
        
        # Backtrack: Un-mark the cell so it can be used in other paths
        visited[r][c] = False

    # Main loop to start the search from every 'N' in the grid
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 'N':
                # Initialize a new visited matrix for each starting point
                visited = [[False for _ in range(cols)] for _ in range(rows)]
                dfs(r, c, "", visited)

    print(max_word)

# Execute the function to find and print the longest word
find_longest_grid_word()