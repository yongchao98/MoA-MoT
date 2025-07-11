import requests
import sys

class TrieNode:
    """A node in the Trie data structure."""
    def __init__(self):
        self.children = {}
        self.is_end_of_word = False

class Trie:
    """Trie data structure for efficient prefix-based searching."""
    def __init__(self):
        self.root = TrieNode()

    def insert(self, word):
        node = self.root
        for char in word:
            if char not in node.children:
                node.children[char] = TrieNode()
            node = node.children[char]
        node.is_end_of_word = True

def solve_word_grid():
    """
    Finds the longest word in a grid starting with a specific letter.
    """
    # Step 1: Define the grid
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    grid_rows, grid_cols = 4, 4
    start_char = 'N'

    # Step 2: Get a word list and build the Trie
    try:
        # Using a well-known public dictionary file.
        # This step requires an internet connection.
        url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
        print("Downloading dictionary file... (this may take a moment)", file=sys.stderr)
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an exception for bad status codes
        words = response.text.splitlines()
        print("Dictionary downloaded.", file=sys.stderr)
    except requests.exceptions.RequestException as e:
        print(f"Error: Could not download the dictionary file. Please check your internet connection.", file=sys.stderr)
        print("Using a small fallback dictionary.", file=sys.stderr)
        # A small list of fallback words if the download fails
        words = ["NOPE", "NOPES", "NOSE", "NOTE", "SON", "SPOKE", "SPOKEN"]

    trie = Trie()
    for word in words:
        if len(word) >= 3: # Ignore very short words
            trie.insert(word.upper())

    # Step 3: Global variables to store the longest word found
    longest_word_found = ""
    longest_word_path = []

    def dfs(r, c, visited, current_path, trie_node):
        """
        Recursively search for words using DFS.
        """
        nonlocal longest_word_found, longest_word_path
        
        char = grid[r][c]
        
        # Move to the next node in the trie
        if char not in trie_node.children:
            return
            
        new_trie_node = trie_node.children[char]
        
        # Form the new word and path
        current_word = "".join(p[0] for p in current_path) + char
        new_path = current_path + [(char, (r, c))]
        
        # If we formed a valid word, check if it's the longest
        if new_trie_node.is_end_of_word:
            if len(current_word) > len(longest_word_found):
                longest_word_found = current_word
                longest_word_path = new_path
        
        # Mark current cell as visited for this path
        visited.add((r, c))
        
        # Explore all 8 neighbors
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                
                nr, nc = r + dr, c + dc
                
                if 0 <= nr < grid_rows and 0 <= nc < grid_cols and (nr, nc) not in visited:
                    dfs(nr, nc, visited, new_path, new_trie_node)
        
        # Backtrack: un-mark the cell as visited to explore other paths
        visited.remove((r, c))

    # Step 4: Start DFS from each 'N' on the grid
    print("Searching for the longest word...", file=sys.stderr)
    for r in range(grid_rows):
        for c in range(grid_cols):
            if grid[r][c] == start_char:
                # The root of the trie doesn't represent a character.
                # The DFS starts from the children of the root.
                dfs(r, c, set(), [], trie.root)

    # Step 5: Print the final result
    if longest_word_found:
        equation_parts = [f"{char}({r},{c})" for char, (r, c) in longest_word_path]
        print(" + ".join(equation_parts) + f" = {longest_word_found}")
    else:
        print(f"No word starting with '{start_char}' could be found.")

if __name__ == '__main__':
    solve_word_grid()