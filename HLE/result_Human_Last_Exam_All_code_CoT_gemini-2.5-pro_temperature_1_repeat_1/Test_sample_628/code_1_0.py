def solve_word_grid():
    """
    Finds the longest word in a 4x4 grid starting with 'N'.
    """
    # The 4x4 grid of letters
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]

    # A dictionary of possible English words.
    # In a real-world application, this would be much larger.
    word_list = [
        "NO", "SO", "PEN", "PENS", "LENS", "SPOKE", "SPOOL", "OPEN", "PONY",
        "NOPE", "NOES", "NOSE", "NOTE", "SONY", "SPOKEN", "JOKE", "LONE",
        "POKE", "POSE", "SEND", "LEND"
    ]

    # Use sets for efficient O(1) average time complexity lookups
    word_set = set(word_list)
    # Create a set of all prefixes for efficient search pruning
    prefix_set = {word[:i] for word in word_set for i in range(1, len(word) + 1)}

    R, C = len(grid), len(grid[0])
    found_words = []

    def dfs(r, c, current_word, path):
        """
        A recursive Depth-First Search function to find words.
        
        Args:
            r (int): Current row.
            c (int): Current column.
            current_word (str): The word formed so far.
            path (set): A set of (row, col) tuples of visited cells.
        """
        # Append the new letter to the current word
        current_word += grid[r][c]

        # If the current string is not a prefix of any word, stop this path
        if current_word not in prefix_set:
            return

        # If the current string is a valid word, add it to our results
        if current_word in word_set:
            found_words.append(current_word)

        # Explore all 8 neighbors (horizontally, vertically, and diagonally)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                # Skip the current cell itself
                if dr == 0 and dc == 0:
                    continue

                nr, nc = r + dr, c + dc

                # Check if the neighbor is within grid bounds and not already visited
                if 0 <= nr < R and 0 <= nc < C and (nr, nc) not in path:
                    new_path = path.copy()
                    new_path.add((nr, nc))
                    dfs(nr, nc, current_word, new_path)

    # Find all starting positions for the letter 'N'
    start_positions = []
    for r in range(R):
        for c in range(C):
            if grid[r][c] == 'N':
                start_positions.append((r, c))

    # Run DFS from each starting 'N'
    for r_start, c_start in start_positions:
        dfs(r_start, c_start, "", {(r_start, c_start)})
    
    # Determine the longest word found
    if not found_words:
        longest_word = "No word starting with 'N' could be found."
    else:
        longest_word = max(found_words, key=len)

    print(f"The longest word found is: {longest_word}")

# Execute the solver function
solve_word_grid()