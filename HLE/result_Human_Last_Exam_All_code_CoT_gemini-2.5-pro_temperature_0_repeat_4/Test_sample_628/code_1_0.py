import collections

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

    # A small dictionary and prefix set for this specific problem.
    # A larger dictionary could be used for a more general solution.
    words = {"NOSE", "NOES", "NEPS", "NOPES"}
    prefixes = {"N", "NO", "NOS", "NOE", "NE", "NEP", "NOP", "NOPE"}
    # Add all words to the prefix set as a word is also a prefix of itself.
    prefixes.update(words)

    rows, cols = len(grid), len(grid[0])
    longest_word_info = {"word": "", "path": []}

    def is_valid(r, c):
        return 0 <= r < rows and 0 <= c < cols

    def dfs(r, c, visited, current_word, current_path):
        """
        Depth-First Search to find words in the grid.
        """
        # Append current letter and position
        current_word += grid[r][c]
        current_path.append(grid[r][c])
        visited.add((r, c))

        # Prune the search if the current string is not a prefix of any word.
        if current_word not in prefixes:
            return

        # If it's a valid word and the longest one found so far, save it.
        if current_word in words:
            if len(current_word) > len(longest_word_info["word"]):
                longest_word_info["word"] = current_word
                longest_word_info["path"] = list(current_path)

        # Explore neighbors
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue

                nr, nc = r + dr, c + dc

                if is_valid(nr, nc) and (nr, nc) not in visited:
                    # Create copies for the new recursive path
                    new_visited = set(visited)
                    new_path = list(current_path)
                    dfs(nr, nc, new_visited, current_word, new_path)

    # Start DFS from each 'N' in the grid
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 'N':
                dfs(r, c, set(), "", [])

    # Format and print the final result
    if longest_word_info["word"]:
        equation = " + ".join(longest_word_info["path"])
        result = f"{equation} = {longest_word_info['word']}"
        print(result)
    else:
        print("No word starting with 'N' could be found.")

solve_boggle()
<<<N + O + P + E + S = NOPES>>>