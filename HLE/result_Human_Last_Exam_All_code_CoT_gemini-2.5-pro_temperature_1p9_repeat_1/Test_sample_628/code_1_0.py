def solve_boggle():
    """
    Finds the longest word starting with 'N' in a 4x4 grid.
    """
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]

    # A built-in dictionary of words to check against.
    # We only need to consider words that start with 'N'.
    words = {"NO", "NOPE", "NOPES", "NOSEY", "NEST", "NET"}

    # For efficient pruning, create a set of all valid prefixes.
    prefixes = set()
    for word in words:
        for i in range(1, len(word) + 1):
            prefixes.add(word[:i])

    rows, cols = len(grid), len(grid[0])
    longest_word = ""
    longest_word_path = []

    def dfs(r, c, current_path, current_word):
        nonlocal longest_word, longest_word_path

        # Prune search if the current word is not a valid prefix
        if current_word not in prefixes:
            return

        # If it's a valid word, check if it's the longest one found so far
        if current_word in words:
            if len(current_word) > len(longest_word):
                longest_word = current_word
                longest_word_path = current_path

        # Explore all 8 neighbors (including diagonals)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue

                nr, nc = r + dr, c + dc

                # Check if the neighbor is valid
                if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in current_path:
                    dfs(nr, nc, current_path + [(nr, nc)], current_word + grid[nr][nc])

    # Start DFS from every 'N' in the grid
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 'N':
                dfs(r, c, [(r, c)], "N")

    # Format and print the final result
    if longest_word:
        path_str_parts = []
        for r, c in longest_word_path:
            letter = grid[r][c]
            path_str_parts.append(f"{letter}({r}, {c})")
        
        result_str = " + ".join(path_str_parts) + f" = {longest_word}"
        print("Found the longest word starting with 'N':")
        print(result_str)
    else:
        print("No words starting with 'N' could be formed from the dictionary.")

solve_boggle()
<<<NOPES>>>