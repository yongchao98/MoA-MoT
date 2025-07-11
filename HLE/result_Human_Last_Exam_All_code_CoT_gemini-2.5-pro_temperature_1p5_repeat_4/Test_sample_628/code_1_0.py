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
    rows, cols = 4, 4

    # A dictionary of valid English words starting with 'N'.
    # This list is curated to contain potential words from the grid.
    word_list = {
        "NEON", "NEONS", "NEP", "NOPE", "NOPES", "NOSE", "NOSEY", "NONES"
    }

    # Create a set of prefixes for efficient pruning during search.
    prefixes = set()
    for word in word_list:
        for i in range(1, len(word) + 1):
            prefixes.add(word[:i])

    # Variables to store the result
    longest_word_found = ""

    def is_valid(r, c):
        """Checks if a coordinate is within the grid bounds."""
        return 0 <= r < rows and 0 <= c < cols

    def search(path):
        """
        A recursive DFS function to find words.
        'path' is a list of (row, col) tuples representing the current path.
        """
        nonlocal longest_word_found
        
        current_word = "".join(grid[r][c] for r, c in path)

        if current_word not in prefixes:
            return

        if current_word in word_list:
            if len(current_word) > len(longest_word_found):
                longest_word_found = current_word

        # Get the last cell in the path to find its neighbors
        last_r, last_c = path[-1]

        # Explore all 8 neighbors (including diagonals)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                
                nr, nc = last_r + dr, last_c + dc

                if is_valid(nr, nc) and (nr, nc) not in path:
                    # Continue the search with the new cell added to the path
                    search(path + [(nr, nc)])

    # Find all starting positions (cells with 'N')
    start_positions = []
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 'N':
                start_positions.append((r, c))

    # Start the search from each 'N'
    for r_start, c_start in start_positions:
        search([(r_start, c_start)])

    # Print the result in the requested equation format
    if longest_word_found:
        equation = " + ".join(list(longest_word_found))
        print(f"{equation} = {longest_word_found}")
    else:
        print("No word starting with 'N' could be found.")

solve_word_grid()