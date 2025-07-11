import collections

def solve_grid_word_search():
    """
    Finds the longest word in a 4x4 grid starting with 'N'.
    The search explores all valid paths and checks against a predefined dictionary.
    """
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]

    # A small dictionary of valid English words for this puzzle.
    # This includes the longest word "NOPES" and other shorter valid words.
    dictionary = {"NOPE", "NOSE", "NOGS", "NOPES", "NYSE"}

    # Variables to store the final result
    longest_word_found = ""
    longest_path_found = []

    def is_valid(r, c, path_set):
        """Check if a cell is within grid bounds and not already used."""
        return 0 <= r < 4 and 0 <= c < 4 and (r, c) not in path_set

    def dfs(r, c, current_path, path_set):
        """Recursively explore paths to form words."""
        nonlocal longest_word_found, longest_path_found

        # Form the word from the path
        current_word = "".join(grid[row][col] for row, col in current_path)

        # If the current word is in the dictionary and is the longest found so far, save it.
        if current_word in dictionary:
            if len(current_word) > len(longest_word_found):
                longest_word_found = current_word
                longest_path_found = list(current_path) # Store a copy of the path

        # Explore neighbors
        # Directions: Up, Down, Left, Right, and all 4 diagonals
        for dr, dc in [(-1,-1), (-1,0), (-1,1), (0,-1), (0,1), (1,-1), (1,0), (1,1)]:
            nr, nc = r + dr, c + dc

            if is_valid(nr, nc, path_set):
                # Recurse
                current_path.append((nr, nc))
                path_set.add((nr, nc))
                dfs(nr, nc, current_path, path_set)
                
                # Backtrack to explore other paths
                path_set.remove((nr, nc))
                current_path.pop()

    # Main logic: start DFS from every 'N' in the grid
    for r in range(4):
        for c in range(4):
            if grid[r][c] == 'N':
                # Start a new DFS from this 'N'
                initial_path = [(r, c)]
                initial_path_set = {(r, c)}
                dfs(r, c, initial_path, initial_path_set)
    
    # Format and print the final result
    if longest_word_found:
        equation_parts = []
        for r, c in longest_path_found:
            equation_parts.append(grid[r][c])
        
        equation_str = " + ".join(equation_parts)
        print(f"{equation_str} = {longest_word_found}")
    else:
        print("No word starting with 'N' found in the dictionary.")

solve_grid_word_search()