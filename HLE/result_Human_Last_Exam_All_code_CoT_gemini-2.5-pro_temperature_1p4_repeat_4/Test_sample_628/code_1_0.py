def solve_word_grid():
    """
    Finds the longest word in a 4x4 grid starting with a specific letter
    using a depth-first search algorithm.
    """
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    start_char = 'N'
    
    # A self-contained dictionary for this specific problem.
    # This makes the code runnable without external files.
    word_set = {"NOPE", "NOSE", "NOSY", "NEPOS", "NOPES"}
    
    # Create a set of all possible prefixes for efficient searching.
    # A path is abandoned if its string is not in this prefix set.
    prefix_set = set()
    for word in word_set:
        for i in range(1, len(word) + 1):
            prefix_set.add(word[:i])

    found_words = []

    def search(r, c, visited, path):
        """
        Recursively explores the grid to find words.
        
        Args:
            r (int): Current row.
            c (int): Current column.
            visited (set): A set of (row, col) tuples already in the path.
            path (list): A list of characters forming the current word.
        """
        # Add the current cell's character and position to the path
        path.append(grid[r][c])
        visited.add((r, c))
        
        current_word = "".join(path)

        # If the current string is not a valid prefix, stop this path
        if current_word not in prefix_set:
            visited.remove((r, c))
            path.pop()
            return

        # If it's a valid word in our dictionary, save it
        if current_word in word_set:
            found_words.append(list(path)) # Add a copy of the path

        # Explore all 8 neighbors
        for dr in range(-1, 2):
            for dc in range(-1, 2):
                if dr == 0 and dc == 0:
                    continue
                
                nr, nc = r + dr, c + dc

                # Check if the neighbor is within bounds and not visited
                if 0 <= nr < len(grid) and 0 <= nc < len(grid[0]) and (nr, nc) not in visited:
                    search(nr, nc, visited, path)
        
        # Backtrack to explore other paths
        visited.remove((r, c))
        path.pop()

    # Start the search from every cell containing the start_char
    for r in range(len(grid)):
        for c in range(len(grid[0])):
            if grid[r][c] == start_char:
                search(r, c, set(), [])

    # Determine the longest word found
    if not found_words:
        print(f"No words starting with '{start_char}' found in the grid.")
    else:
        longest_path = max(found_words, key=len)
        longest_word = "".join(longest_path)
        
        # Format the output as an equation
        equation = " + ".join(longest_path)
        print(f"{equation} = {longest_word}")

solve_word_grid()
<<<N + O + P + E + S = NOPES>>>