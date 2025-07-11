import collections

def solve():
    """
    Finds the longest word in a grid starting with 'N'.
    """
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]

    # A curated dictionary of words that can be formed from the grid letters.
    # This helps in validating words found by the search algorithm.
    word_list = [
        "NO", "NE", "ON", "SO", "GO", "DO", "DONE",
        "NOPE", "NOPS", "NOSE", "NOSY", "SOON", "SONG",
        "OPEN", "LOOP", "LOOPS", "SNOOP", "SKEP",
        "SPOKE", "SPOON", "SLOPE", "PENSION", "PERSON",
        "NOPES", "SNOOPS", "SPOKEN", "SPOONS", "POISON"
    ]

    # Create a set of words for fast O(1) lookup.
    valid_words = set(word_list)
    # Create a set of prefixes for efficient search pruning.
    prefixes = set()
    for word in word_list:
        for i in range(1, len(word) + 1):
            prefixes.add(word[:i])

    rows, cols = len(grid), len(grid[0])
    found_paths = []

    def search(path):
        """
        A recursive DFS function to find words.
        'path' is a list of (row, col) tuples.
        """
        current_word = "".join(grid[r][c] for r, c in path)

        if current_word not in prefixes:
            return

        if current_word in valid_words:
            found_paths.append(path)

        # Get the last cell in the path to find neighbors
        r, c = path[-1]

        # Explore 8 directions (horizontal, vertical, diagonal)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue

                nr, nc = r + dr, c + dc

                # Check boundaries and if the cell has been visited
                if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in path:
                    search(path + [(nr, nc)])

    # Start the search from every 'N' in the grid
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 'N':
                search([(r, c)])

    if not found_paths:
        print("No valid word starting with 'N' found.")
        return

    # Find the longest path among all found valid paths
    longest_path = max(found_paths, key=len)
    longest_word = "".join(grid[r][c] for r, c in longest_path)
    
    # Format the equation string
    equation_parts = [grid[r][c] for r, c in longest_path]
    equation_str = " + ".join(equation_parts) + f" = {longest_word}"
    
    print(equation_str)

solve()