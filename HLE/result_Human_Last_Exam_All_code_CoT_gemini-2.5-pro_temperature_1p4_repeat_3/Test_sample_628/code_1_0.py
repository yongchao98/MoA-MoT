import collections

def solve_word_grid():
    """
    Finds the longest word in a grid starting with 'N' using DFS and a Trie.
    """
    # The 4x4 grid of letters
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]

    # A dictionary of words to check against.
    # A Trie will be built from these words for efficient prefix checking.
    # All words are in uppercase to match the grid.
    WORDS = {
        "NO", "NOPE", "NOPES", "NOSE", "NOTE", "NOTES", "NYES", "POSE",
        "KENO", "OPEN", "PENS", "SON", "SOY", "SPY", "YEN", "YES", "EON",
        "EONS", "GO", "PONY"
    }

    # Define the TrieNode for our prefix tree
    class TrieNode:
        def __init__(self):
            self.children = collections.defaultdict(TrieNode)
            self.is_word = False

    # Create the root of the Trie
    trie_root = TrieNode()

    # Populate the Trie with our dictionary of words
    for word in WORDS:
        node = trie_root
        for char in word:
            node = node.children[char]
        node.is_word = True

    # This set will store all valid words found in the grid
    found_words = set()

    # The Depth-First Search function to traverse the grid
    def dfs(r, c, node, current_word, visited):
        """
        Recursively explores the grid to find words.
        :param r: current row
        :param c: current column
        :param node: current node in the Trie
        :param current_word: the word formed so far
        :param visited: a set of (row, col) tuples of visited cells for the current path
        """
        # If the current path of letters forms a valid word, add it to our set
        if node.is_word:
            found_words.add(current_word)

        # Explore all 8 neighbors (horizontal, vertical, diagonal)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue  # Skip the cell itself

                nr, nc = r + dr, c + dc

                # Check if the neighbor is within grid boundaries and hasn't been visited yet
                if 0 <= nr < 4 and 0 <= nc < 4 and (nr, nc) not in visited:
                    char = grid[nr][nc]
                    # If the next character continues a valid prefix in the Trie
                    if char in node.children:
                        visited.add((nr, nc))
                        dfs(nr, nc, node.children[char], current_word + char, visited)
                        # Backtrack: un-visit the cell to allow it to be used in other paths
                        visited.remove((nr, nc))

    # --- Main Search Logic ---
    # Iterate through the grid to find starting letters
    for r in range(4):
        for c in range(4):
            start_char = grid[r][c]
            # We only start the search from the letter 'N'
            if start_char == 'N':
                # Check if 'N' can be the start of any word in our dictionary
                if start_char in trie_root.children:
                    # Start the DFS from this cell
                    visited = set([(r, c)])
                    dfs(r, c, trie_root.children[start_char], start_char, visited)

    # After the search is complete, find the longest word
    if not found_words:
        longest_word = "No valid word found starting with 'N'."
    else:
        # Find the word with the maximum length from the set of found words
        longest_word = max(found_words, key=len)

    # --- Final Output ---
    print("Grid:")
    for row in grid:
        print(f"  {' '.join(row)}")
    print("\nSearching for the longest word starting with 'N'...")
    
    if found_words:
        print(f"Found {len(found_words)} valid word(s): {sorted(list(found_words))}")
        print(f"\nThe longest word is: {longest_word}")
    else:
        print("\nNo valid words starting with 'N' could be formed.")

# Run the solver
solve_word_grid()
<<<NOPES>>>