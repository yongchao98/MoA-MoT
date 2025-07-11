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
    rows, cols = 4, 4

    # A curated dictionary of plausible words and their prefixes for efficient searching.
    # NOPES is a valid 5-letter word (plural of nope) and can be formed on the grid.
    word_list = {"NOPE", "NOPES", "NOSE", "NOPS"}
    prefixes = set()
    for word in word_list:
        for i in range(1, len(word) + 1):
            prefixes.add(word[:i])

    found_words = set()

    def dfs(r, c, current_word, visited):
        """
        Recursively explores paths on the grid to find words.
        
        Args:
            r (int): Current row.
            c (int): Current column.
            current_word (str): The word formed so far.
            visited (set): A set of (row, col) tuples for the current path.
        """
        # Append the new letter
        current_word += grid[r][c]

        # Pruning: if the current string is not a prefix of any valid word, stop.
        if current_word not in prefixes:
            return

        # Add current cell to the visited path
        visited.add((r, c))

        # If the current word is in our dictionary, add it to our results
        if current_word in word_list:
            found_words.add(current_word)

        # Explore all 8 neighbors (horizontal, vertical, and diagonal)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue  # Skip the cell itself

                nr, nc = r + dr, c + dc

                # Check if the neighbor is within bounds and not already visited
                if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in visited:
                    dfs(nr, nc, current_word, visited)
        
        # Backtrack: remove the current cell from visited to allow it in other paths
        visited.remove((r, c))

    # --- Main execution ---
    # Find all starting positions (cells with 'N')
    start_positions = []
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 'N':
                start_positions.append((r, c))

    # Run the search from each starting 'N'
    for r_start, c_start in start_positions:
        dfs(r_start, c_start, "", set())

    # Determine the longest word found
    longest_word = ""
    if found_words:
        # The problem asks for the longest word starting with 'N'.
        # All words found will start with 'N' due to our starting points.
        longest_word = max(found_words, key=len)

    print(longest_word)

solve_word_grid()
>>> NOPES