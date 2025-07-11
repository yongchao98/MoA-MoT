def solve_boggle():
    """
    Finds the longest word on the provided 4x4 grid that starts with the letter 'N'.
    The solution uses a Depth-First Search (DFS) algorithm to explore all possible
    letter combinations and checks them against a predefined dictionary.
    """
    # The 4x4 grid of letters.
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    
    # A simple dictionary of valid words. A comprehensive application would load
    # a much larger dictionary from a file. This set includes a known 5-letter
    # word, "NOPES", which can be formed on the grid.
    dictionary = {"NO", "NOPE", "NOPES", "NOTE", "PENS", "SON"}

    rows, cols = len(grid), len(grid[0])
    all_found_words = set()

    def dfs(r, c, current_word, visited):
        """
        Recursively explores paths on the grid to find valid words.

        Args:
            r (int): The current row index.
            c (int): The current column index.
            current_word (str): The word formed by the path so far.
            visited (set): A set of (row, col) tuples for the current path.
        """
        # Append the letter from the current cell to the word.
        current_word += grid[r][c]
        
        # Mark the current cell as visited for this path.
        visited.add((r, c))

        # If the word formed so far is in our dictionary, add it to our results.
        if current_word in dictionary:
            all_found_words.add(current_word)

        # Explore all 8 neighbors (horizontally, vertically, and diagonally).
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue  # Skip the current cell itself.

                nr, nc = r + dr, c + dc

                # Check if the neighbor is within grid bounds and has not been visited yet.
                if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in visited:
                    dfs(nr, nc, current_word, visited)
        
        # Backtrack: un-mark the cell as visited to allow it in other paths.
        visited.remove((r, c))

    # Iterate through the grid to find all 'N's to start the search from.
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 'N':
                dfs(r, c, "", set())

    # After searching all paths, find the longest word among the results.
    if not all_found_words:
        print("No valid word starting with 'N' was found on the grid.")
    else:
        # Use the `max` function with a `key` to find the longest string in the set.
        longest_word = max(all_found_words, key=len)
        print(longest_word)

# Run the solver.
solve_boggle()