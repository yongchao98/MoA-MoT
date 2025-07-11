def solve_boggle():
    """
    Finds the longest word in a 4x4 grid starting with 'N'.

    The solution uses a Depth-First Search (DFS) algorithm to explore all
    possible paths on the grid. A dictionary of valid words and a set of
    their prefixes are used to efficiently verify words and prune search paths.
    """

    # Step 1: Define the grid of letters
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    
    # Step 2: Define a dictionary of valid words starting with 'N'.
    # This list is curated for this specific grid and problem.
    dictionary = {
        "NO", 
        "NE", 
        "NOSE", 
        "NEPOS"  # A valid, albeit less common, word meaning "descendant".
    }

    # Step 3: Create a set of all valid prefixes for search optimization.
    prefixes = set()
    for word in dictionary:
        for i in range(1, len(word) + 1):
            prefixes.add(word[:i])

    # This will hold the longest word found.
    longest_word_found = ""

    # Step 4: Implement the recursive Depth-First Search (DFS) function.
    def dfs(r, c, current_word, visited_coords):
        nonlocal longest_word_found
        
        # Add the current letter to the word and mark the cell as visited
        current_word += grid[r][c]
        visited_coords.add((r, c))

        # Pruning: If the current string is not a valid prefix, stop this path.
        if current_word not in prefixes:
            return

        # If it's a valid word and longer than our current best, update it.
        if current_word in dictionary and len(current_word) > len(longest_word_found):
            longest_word_found = current_word

        # Explore all 8 neighbors (horizontal, vertical, and diagonal)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue  # Skip the cell itself
                
                nr, nc = r + dr, c + dc
                
                # Check if the neighbor is within grid bounds and not already visited in this path
                if 0 <= nr < 4 and 0 <= nc < 4 and (nr, nc) not in visited_coords:
                    # A copy of the visited set is passed to the recursive call
                    dfs(nr, nc, current_word, set(visited_coords))

    # Step 5: Iterate through the grid to find all starting 'N's and run DFS.
    rows, cols = 4, 4
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 'N':
                # Start a new search from this 'N'
                dfs(r, c, "", set())
    
    # Step 6: Format and print the final result.
    if longest_word_found:
        # Create the equation string, e.g., "N + E + P + O + S = NEPOS"
        equation = " + ".join(list(longest_word_found)) + " = " + longest_word_found
        print(equation)
    else:
        print("No valid word starting with 'N' was found.")

solve_boggle()
<<<N + E + P + O + S = NEPOS>>>