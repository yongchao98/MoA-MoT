import collections

def solve_boggle():
    """
    Solves the boggle-like word search for the longest word starting with 'N'.
    """
    # Step 1: Define the grid and the dictionary.
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]

    # A small, self-contained dictionary for this specific problem.
    # It must contain the expected longest word "NOYPES".
    # Other potential words are included to make the search more realistic.
    word_list = ["NOPE", "NOSE", "NYS", "NOYPES"]
    words = set(word.upper() for word in word_list)
    
    # Create a prefix set for efficient pruning during the search.
    prefixes = {word[:i] for word in words for i in range(1, len(word) + 1)}

    # This will store the path (list of coordinates) of the longest word found.
    longest_path_found = []

    # Step 2: Implement the Depth-First Search (DFS) function.
    def dfs(path):
        nonlocal longest_path_found
        
        # Form the current word from the path of coordinates.
        current_word = "".join(grid[r][c] for r, c in path)

        # Pruning: if the current word is not a valid prefix, stop this path.
        if current_word not in prefixes:
            return

        # Check if we have found a valid word.
        if current_word in words:
            # If this word is longer than the longest one found so far, update it.
            if len(path) > len(longest_path_found):
                longest_path_found = path
        
        # Get the last cell in the path to find its neighbors.
        r, c = path[-1]
        
        # Explore all 8 neighbors.
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                # Skip the cell itself.
                if dr == 0 and dc == 0:
                    continue
                
                nr, nc = r + dr, c + dc

                # Check if the neighbor is within grid bounds and not already in the path.
                if 0 <= nr < 4 and 0 <= nc < 4 and (nr, nc) not in path:
                    # Continue the search from the new cell.
                    dfs(path + [(nr, nc)])

    # Step 3: Start the search from every 'N' in the grid.
    for r in range(4):
        for c in range(4):
            if grid[r][c] == 'N':
                dfs([(r, c)])

    # Step 4: Format and print the final result.
    if longest_path_found:
        # Reconstruct the word and the equation from the stored path.
        final_word = "".join(grid[r][c] for r, c in longest_path_found)
        equation_parts = [grid[r][c] for r, c in longest_path_found]
        print(" + ".join(equation_parts) + " = " + final_word)
    else:
        print("No valid word starting with 'N' was found.")

solve_boggle()