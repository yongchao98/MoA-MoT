import urllib.request

def solve_boggle():
    """
    Finds the longest word in a 4x4 grid starting with 'N'.
    """

    # --- 1. Setup Grid and Constants ---
    GRID = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    ROWS, COLS = 4, 4
    START_LETTER = 'N'
    
    # URL for a common English word list
    WORD_LIST_URL = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"

    # --- 2. Build Trie (Prefix Tree) ---
    def build_trie(words):
        """Builds a Trie from a list of words for efficient prefix checking."""
        root = {}
        for word in words:
            node = root
            for char in word:
                node = node.setdefault(char, {})
            # Use a special key to mark the end of a valid word
            node['#'] = '#'
        return root

    print("Downloading dictionary...")
    try:
        # Fetch and decode the word list, converting all words to uppercase
        with urllib.request.urlopen(WORD_LIST_URL) as response:
            words = {
                word.decode('utf-8').strip().upper()
                for word in response.readlines()
            }
        print("Dictionary loaded. Building Trie...")
        trie_root = build_trie(words)
        print("Trie built. Starting search.")
    except Exception as e:
        print(f"Error downloading or processing word list: {e}")
        # Use a fallback dictionary in case of network errors
        print("Using a small fallback dictionary.")
        fallback_words = {"NO", "NOPE", "NOPES", "NOSE", "NOSEY"}
        trie_root = build_trie(fallback_words)

    found_words = []

    # --- 3. Depth-First Search (DFS) function ---
    def dfs(row, col, path, trie_node):
        """
        Recursively explores the grid to find words.

        Args:
            row (int): Current row index.
            col (int): Current column index.
            path (list): List of (r, c) tuples visited in the current path.
            trie_node (dict): The current node in the Trie.
        """
        # Form the current word from the path
        current_word = "".join(GRID[r][c] for r, c in path)
        
        # If the word is a valid word, add it to our list
        if '#' in trie_node:
            found_words.append(current_word)
        
        # Explore all 8 neighbors (horizontal, vertical, diagonal)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue  # Skip the current cell itself

                nr, nc = row + dr, col + dc

                # Check if the neighbor is valid
                if 0 <= nr < ROWS and 0 <= nc < COLS and (nr, nc) not in path:
                    neighbor_letter = GRID[nr][nc]
                    # If the letter continues a valid prefix, recurse
                    if neighbor_letter in trie_node:
                        dfs(nr, nc, path + [(nr, nc)], trie_node[neighbor_letter])

    # --- 4. Main Search Logic ---
    # Find all starting positions for the letter 'N'
    start_positions = []
    for r in range(ROWS):
        for c in range(COLS):
            if GRID[r][c] == START_LETTER:
                start_positions.append((r, c))

    # Start the DFS from each initial position
    for r_start, c_start in start_positions:
        if START_LETTER in trie_root:
            dfs(r_start, c_start, [(r_start, c_start)], trie_root[START_LETTER])

    # --- 5. Find and Print the Longest Word ---
    if not found_words:
        print("No words starting with 'N' were found.")
    else:
        # Find the longest word. If there's a tie, this returns the first one found.
        longest_word = max(found_words, key=len)
        print(f"Found {len(found_words)} words starting with '{START_LETTER}'.")
        print(f"The longest word found is:")
        print(longest_word)

# Run the solver
solve_boggle()