import sys

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

    # A fallback dictionary in case /usr/share/dict/words is not available
    FALLBACK_WORDS = [
        "NO", "SO", "SON", "NE", "ON", "OPEN", "OPENS", "OPES", "SPY",
        "NOPE", "NOES", "NOSE", "NOSY", "POKE", "POKES", "PONY", "POSE",
        "SONY", "SNOOP", "SNOOPY", "SPOKE", "SPOKEN", "NEON", "OPES",
        "PENS", "PENSY", "SOON", "SPOON", "NOPES"
    ]

    def load_dictionary():
        """
        Loads words from a dictionary file and creates a prefix set for optimization.
        Returns a set of valid words and a set of valid prefixes.
        """
        words = set()
        try:
            with open('/usr/share/dict/words', 'r') as f:
                for line in f:
                    words.add(line.strip().upper())
        except FileNotFoundError:
            print("'/usr/share/dict/words' not found. Using a small fallback dictionary.", file=sys.stderr)
            words = {word.upper() for word in FALLBACK_WORDS}

        prefixes = set()
        for word in words:
            for i in range(1, len(word) + 1):
                prefixes.add(word[:i])
        return words, prefixes

    words, prefixes = load_dictionary()
    rows, cols = len(grid), len(grid[0])
    found_words = set()

    def dfs(r, c, path, current_word):
        """
        Performs a depth-first search to find words on the grid.
        
        Args:
            r (int): Current row index.
            c (int): Current column index.
            path (set): A set of (row, col) tuples representing visited cells.
            current_word (str): The word formed so far.
        """
        # Pruning step: if the current word is not a prefix of any valid word, stop this path.
        if current_word not in prefixes:
            return

        # If the current word is in our dictionary, add it to our results.
        if current_word in words:
            found_words.add(current_word)

        # Explore all 8 neighbors (horizontal, vertical, and diagonal)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue  # Skip the cell itself

                nr, nc = r + dr, c + dc

                # Check if the neighbor is within grid bounds and hasn't been visited in the current path
                if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in path:
                    new_path = path.copy()
                    new_path.add((nr, nc))
                    dfs(nr, nc, new_path, current_word + grid[nr][nc])

    # Start the DFS from every cell containing 'N'
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 'N':
                dfs(r, c, {(r, c)}, 'N')

    # Find and print the longest word among all found words
    if not found_words:
        print("No word starting with 'N' could be formed.")
    else:
        longest_word = max(found_words, key=len)
        print(longest_word)

solve_boggle()