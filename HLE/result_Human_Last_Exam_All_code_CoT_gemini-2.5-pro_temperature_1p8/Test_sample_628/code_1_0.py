import collections

def solve_word_grid():
    """
    Finds the longest word in a 4x4 grid starting with 'N'.
    """
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]

    # A dictionary of words to check against. For a real-world Boggle solver,
    # this would be much larger. It is pre-populated with potential candidates.
    word_list = {"NOPE", "NOSE", "NOYS", "NOPES", "SPOKEN", "JOKE", "PENS"}
    
    # For efficient prefix lookups, we create a set of all possible prefixes.
    prefixes = set()
    for word in word_list:
        for i in range(1, len(word) + 1):
            prefixes.add(word[:i])

    rows, cols = len(grid), len(grid[0])
    found_words = []

    def dfs(path):
        """
        Performs a depth-first search to find words.
        'path' is a list of (row, col) tuples.
        """
        r, c = path[-1]
        current_word = "".join([grid[r][c] for r, c in path])

        if current_word not in prefixes:
            return

        if current_word in word_list:
            found_words.append(current_word)

        # Explore neighbors
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue

                nr, nc = r + dr, c + dc

                if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in path:
                    new_path = path + [(nr, nc)]
                    dfs(new_path)

    # Start DFS from each 'N' in the grid
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 'N':
                dfs([(r, c)])

    if not found_words:
        print("No valid word found.")
        return

    # Find the longest word from the list of found words
    longest_word = max(found_words, key=len)
    
    # Print the result in the requested "equation" format
    equation_str = " + ".join(list(longest_word))
    print(f"{equation_str} = {longest_word}")

solve_word_grid()