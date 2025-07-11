def solve_word_grid():
    """
    Finds the longest word starting with 'N' in a 4x4 grid.
    """
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    
    # A small dictionary for this problem. A real application would use a large word list.
    # Words are included that could plausibly be formed on the grid.
    dictionary = {"NOPE", "NOSE", "POKE", "SPOKE", "LENS", "OPEN", "POESY"}
    
    rows, cols = 4, 4
    longest_word_info = {"word": "", "path": []}

    def is_valid(r, c, visited):
        """Check if a cell is within grid bounds and not visited."""
        return 0 <= r < rows and 0 <= c < cols and not visited[r][c]

    def get_char(r, c):
        """Get character from grid coordinates."""
        return grid[r][c]

    def dfs(r, c, path, visited):
        """Recursive DFS to find words."""
        nonlocal longest_word_info
        
        path.append((r, c))
        visited[r][c] = True
        
        current_word = "".join([get_char(pr, pc) for pr, pc in path])
        
        # Check if the current word is in the dictionary and is the longest found so far.
        if current_word in dictionary:
            if len(current_word) > len(longest_word_info["word"]):
                longest_word_info["word"] = current_word
                longest_word_info["path"] = list(path) # Store a copy of the path

        # Explore all 8 neighbors
        for dr in range(-1, 2):
            for dc in range(-1, 2):
                if dr == 0 and dc == 0:
                    continue
                
                nr, nc = r + dr, c + dc
                
                if is_valid(nr, nc, visited):
                    # A more optimal solution would use a Trie (prefix tree) to prune
                    # searches where the current_word is not a valid prefix.
                    # For this small grid, we explore all paths up to a certain depth.
                    # We can add a simple length check to stop excessively long paths.
                    if len(path) < 8: # Arbitrary limit to prevent infinite loops on faulty logic
                        dfs(nr, nc, path, visited)

        # Backtrack: un-visit the cell and remove it from the path
        visited[r][c] = False
        path.pop()

    # Main logic starts here
    start_char = 'N'
    for r_start in range(rows):
        for c_start in range(cols):
            if get_char(r_start, c_start) == start_char:
                # Initialize visited grid for each new starting 'N'
                visited = [[False for _ in range(cols)] for _ in range(rows)]
                dfs(r_start, c_start, [], visited)

    # Print the final result
    if longest_word_info["word"]:
        equation_letters = [get_char(r, c) for r, c in longest_word_info["path"]]
        print(" + ".join(equation_letters) + " = " + longest_word_info["word"])
    else:
        print(f"No word starting with '{start_char}' found in the dictionary.")

solve_word_grid()
<<<N + O + P + E = NOPE>>>