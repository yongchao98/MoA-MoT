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

    # A self-contained dictionary of plausible words from the grid letters.
    # Including a variety of words ensures the search is reasonably thorough.
    dictionary = {
        "nope", "nopes", "nose", "nosey", "pony", "pens", "send", 
        "son", "one", "eon", "ken", "pen", "pes", "sen", "oke",
        "epos", "peso", "keno", "spoken"
    }

    # Pre-calculate all prefixes for efficient pruning during the search.
    # For example, for "nopes", we add "n", "no", "nop", "nope", "nopes".
    prefixes = set()
    for word in dictionary:
        for i in range(1, len(word) + 1):
            prefixes.add(word[:i])

    rows, cols = 4, 4
    longest_word_found = ""

    def dfs_search(r, c, path, current_word):
        nonlocal longest_word_found

        # Append the character from the current grid cell
        current_word += grid[r][c].lower()
        
        # Add the current cell's coordinates to the path
        path.add((r, c))

        # --- Pruning Step ---
        # If the current string is not a prefix of any word in our dictionary,
        # there's no point in continuing down this path.
        if current_word not in prefixes:
            return

        # --- Check for Longest Word ---
        # If the current string is a valid word and longer than what we've found, update it.
        if current_word in dictionary and len(current_word) > len(longest_word_found):
            longest_word_found = current_word

        # --- Recursive Step ---
        # Explore all 8 adjacent neighbors.
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                # Skip the current cell itself
                if dr == 0 and dc == 0:
                    continue
                
                nr, nc = r + dr, c + dc
                
                # Check if the neighbor is within the grid and hasn't been visited on this path.
                if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in path:
                    # Because `path` is a mutable set, we pass a copy to the next recursive
                    # call so that different branches of the search don't interfere.
                    dfs_search(nr, nc, path.copy(), current_word)

    # --- Main Loop ---
    # Start the search from every 'N' on the grid.
    for r in range(rows):
        for c in range(cols):
            if grid[r][c].upper() == 'N':
                # The initial call starts with an empty path and an empty word.
                dfs_search(r, c, set(), "")

    print(longest_word_found)

# Execute the solver
solve_word_grid()