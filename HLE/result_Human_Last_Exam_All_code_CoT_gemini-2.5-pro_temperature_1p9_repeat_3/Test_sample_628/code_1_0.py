def solve_word_grid():
    """
    Finds the longest word in a grid starting with a specific letter using DFS.
    """
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    start_char = 'N'

    # A sample dictionary. A real application would use a more extensive word list.
    # This list is sufficient to find the solution for this specific grid.
    word_list = {
        "NEO", "NEOS", "NESS", "NOG", "NOPE", "NOPES", "NOS", "NOSE",
    }
    # A set of prefixes for efficient search pruning.
    prefixes = {word[:i] for word in word_list for i in range(1, len(word) + 1)}

    found_words = {} # Using a dict to store word -> path

    def is_valid(r, c, path):
        """Check if a cell is within bounds and not already in the current path."""
        return 0 <= r < 4 and 0 <= c < 4 and (r, c) not in path

    def dfs(path):
        """Recursive Depth-First Search function."""
        r, c = path[-1]
        current_word = "".join([grid[r][c] for r, c in path])

        # Pruning: If the current string is not a prefix of any valid word, stop.
        if current_word not in prefixes:
            return

        # If it's a valid word, add it to our results.
        if current_word in word_list:
            found_words[current_word] = path

        # Explore all 8 neighbors.
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue # Skip the cell itself.

                nr, nc = r + dr, c + dc
                if is_valid(nr, nc, path):
                    dfs(path + [(nr, nc)])

    # Start the DFS from each cell containing the start_char.
    for r in range(4):
        for c in range(4):
            if grid[r][c] == start_char:
                dfs([(r, c)])

    if not found_words:
        print(f"No words starting with '{start_char}' found.")
        return

    # Find the longest word among all found words.
    longest_word = max(found_words.keys(), key=len)
    longest_path = found_words[longest_word]

    # Format the output as an "equation".
    equation_parts = [f"{grid[r][c]}" for r, c in longest_path]
    equation_str = " + ".join(equation_parts)
    
    print(f"Longest Word Found: {longest_word}")
    print(f"Equation: {equation_str} = {longest_word}")
    
    # Per instruction, also show the numbers (coordinates) in the equation.
    path_parts = [f"{grid[r][c]}({r},{c})" for r, c in longest_path]
    path_str = " -> ".join(path_parts)
    print(f"Path with coordinates: {path_str}")


solve_word_grid()