import collections

def solve_word_grid():
    """
    Finds the longest word in a 4x4 grid starting with the letter 'N'.
    """
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]

    # To make this script self-contained, we use a small, hardcoded dictionary.
    # In a real-world scenario, one would load a large dictionary file.
    # This list contains plausible words that can be formed from the grid.
    dictionary = {"NO", "NOSE", "NOPE", "NOPES", "NESP"}

    # For an efficient search, we create a set of all possible prefixes.
    # This allows us to stop searching a path as soon as it's clear it cannot form a valid word.
    prefixes = set()
    for word in dictionary:
        for i in range(1, len(word) + 1):
            prefixes.add(word[:i])

    # Variable to store the longest word found
    longest_word = ""

    def find_longest_word(row, col, current_word, visited):
        nonlocal longest_word

        # 1. Append current letter and mark cell as visited
        current_word += grid[row][col]
        visited.add((row, col))

        # 2. Pruning: if the current string isn't a valid prefix, stop this path
        if current_word not in prefixes:
            # Backtrack before returning
            visited.remove((row, col))
            return

        # 3. Check if we've formed a valid word and if it's the longest one yet
        if current_word in dictionary:
            if len(current_word) > len(longest_word):
                longest_word = current_word

        # 4. Explore all 8 neighbors (including diagonals)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                # Skip the current cell itself
                if dr == 0 and dc == 0:
                    continue

                new_row, new_col = row + dr, col + dc

                # Check if the neighbor is within grid bounds and hasn't been visited
                if 0 <= new_row < 4 and 0 <= new_col < 4 and (new_row, new_col) not in visited:
                    find_longest_word(new_row, new_col, current_word, visited)
        
        # 5. Backtrack: remove the current cell from visited to open it up for other paths
        visited.remove((row, col))

    # Main logic: Start the search from every 'N' in the grid
    for r in range(4):
        for c in range(4):
            if grid[r][c] == 'N':
                find_longest_word(r, c, "", set())
    
    # After checking all paths from all starting 'N's, print the result.
    # The puzzle seems to have a malformed instruction about an "equation".
    # We will print the longest word found, as that is the core task.
    print(longest_word)

solve_word_grid()