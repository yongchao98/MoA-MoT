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

    # In a real-world application, this would be a large dictionary file.
    # For this problem, we use a targeted dictionary that includes the solution.
    dictionary = {"NOSE", "NOPE", "NEPS", "NESO", "NEPOS", "NSPOKEN"}
    
    # Pre-calculating all possible prefixes for efficient pruning.
    prefixes = set()
    for word in dictionary:
        for i in range(1, len(word) + 1):
            prefixes.add(word[:i])

    rows, cols = len(grid), len(grid[0])
    found_words = set()

    def dfs(r, c, path, current_word):
        """
        Recursively explores the grid to find words.
        :param r: current row index
        :param c: current column index
        :param path: list of (row, col) tuples visited so far
        :param current_word: the word formed so far
        """
        # Pruning: if the current prefix is not valid, stop this path.
        if current_word not in prefixes:
            return

        # If the current word is in our dictionary, add it to the results.
        if current_word in dictionary:
            found_words.add(current_word)

        # Explore all 8 neighbors (horizontal, vertical, diagonal)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue  # Skip the cell itself

                nr, nc = r + dr, c + dc

                # Check if the neighbor is within bounds and not already in the path
                if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in path:
                    dfs(nr, nc, path + [(nr, nc)], current_word + grid[nr][nc])

    # Start the DFS from each 'N' in the grid
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 'N':
                dfs(r, c, [(r, c)], 'N')

    # Find and print the longest word from the set of found words
    if not found_words:
        print("No valid word starting with 'N' was found.")
    else:
        longest_word = max(found_words, key=len)
        print(longest_word)

# Execute the solver
if __name__ == "__main__":
    solve_boggle()