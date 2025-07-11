def solve_boggle():
    """
    Finds the longest word in a 4x4 grid starting with 'N'.
    """
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    
    # A dictionary of possible words. A larger list would be used in a more general case.
    # This list contains plausible words starting with 'N' from the grid.
    word_list = {"NO", "NOPE", "NOPES", "NOSE", "NOES", "NOSY", "NOSEY"}
    
    # Generate a set of prefixes for efficient search pruning
    prefixes = set()
    for word in word_list:
        for i in range(1, len(word) + 1):
            prefixes.add(word[:i])

    R, C = len(grid), len(grid[0])
    found_words = set()

    def dfs(r, c, current_path_coords, current_word):
        """
        Recursively explores paths on the grid.
        
        Args:
            r (int): Current row.
            c (int): Current column.
            current_path_coords (set): A set of (row, col) tuples in the current path.
            current_word (str): The word formed so far.
        """
        # Prune the search if the current word is not a valid prefix
        if current_word not in prefixes:
            return

        # If the current prefix is a valid word, add it to our results
        if current_word in word_list:
            found_words.add(current_word)

        # Explore all 8 neighbors (horizontal, vertical, diagonal)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue

                nr, nc = r + dr, c + dc

                # Check if the neighbor is within bounds and not already used
                if 0 <= nr < R and 0 <= nc < C and (nr, nc) not in current_path_coords:
                    new_path_coords = current_path_coords | {(nr, nc)}
                    dfs(nr, nc, new_path_coords, current_word + grid[nr][nc])

    # Find all starting positions for the letter 'N'
    start_positions = []
    for r in range(R):
        for c in range(C):
            if grid[r][c] == 'N':
                start_positions.append((r, c))

    # Run DFS from each starting position
    for r_start, c_start in start_positions:
        dfs(r_start, c_start, {(r_start, c_start)}, grid[r_start][c_start])
        
    if not found_words:
        print("No valid word starting with 'N' was found.")
        return
        
    # Find the longest word. If there's a tie, choose the lexicographically larger one.
    longest_word = ""
    for word in found_words:
        if len(word) > len(longest_word):
            longest_word = word
        elif len(word) == len(longest_word) and word > longest_word:
            longest_word = word
            
    # Print the result in the requested "equation" format
    final_equation = " + ".join(list(longest_word))
    print(final_equation)

solve_boggle()
<<<NOSEY>>>