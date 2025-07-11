import urllib.request
import urllib.error
import sys

class TrieNode:
    """A node in the Trie data structure."""
    def __init__(self):
        self.children = {}
        self.is_end_of_word = False

def build_trie(words):
    """Builds a Trie from a set of words."""
    root = TrieNode()
    for word in words:
        node = root
        for char in word:
            if char not in node.children:
                node.children[char] = TrieNode()
            node = node.children[char]
        node.is_end_of_word = True
    return root

def solve_boggle():
    """
    Finds the longest word in a grid starting with a specific letter.
    """
    # 1. Setup the grid and starting letter
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    start_letter = 'N'
    rows, cols = len(grid), len(grid[0])

    # 2. Get the word list and build the Trie
    words = set()
    try:
        # A common, comprehensive word list
        url = "https://raw.githubusercontent.com/dwyl/english-words/master/words_alpha.txt"
        print("Downloading word list...", file=sys.stderr)
        with urllib.request.urlopen(url) as response:
            # Decode, convert to uppercase, and strip whitespace
            text = response.read().decode('utf-8')
            # Filter for words of length 3 or more
            words = {word.upper() for word in text.splitlines() if len(word) >= 3}
        print("Word list downloaded and processed.", file=sys.stderr)
    except (urllib.error.URLError, ConnectionError) as e:
        print(f"Error: Could not download word list ({e}).", file=sys.stderr)
        print("Using a small fallback dictionary. Results may be limited.", file=sys.stderr)
        # A small fallback dictionary including potential answers
        words = {"NOSE", "SOP", "SON", "PENS", "PENT", "POSE", "OPEN"}
        
    trie_root = build_trie(words)
    
    # 3. DFS implementation
    found_words = set()

    def find_words_recursive(r, c, current_node, path):
        """Recursively search for words from cell (r, c)."""
        letter = grid[r][c]

        # Prune search if the current letter does not form a valid prefix
        if letter not in current_node.children:
            return

        # Move to the next node in the Trie and update the path
        path.append((r, c))
        next_node = current_node.children[letter]

        # If a valid word is found, add it to our set
        if next_node.is_end_of_word:
            word = "".join(grid[i][j] for i, j in path)
            found_words.add(word)

        # Explore all 8 neighbors
        for dr in range(-1, 2):
            for dc in range(-1, 2):
                if dr == 0 and dc == 0:
                    continue # Skip the cell itself

                nr, nc = r + dr, c + dc
                
                # Check for valid and unvisited neighbor
                if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in path:
                    find_words_recursive(nr, nc, next_node, path)
        
        # Backtrack: remove the current cell from the path to explore other branches
        path.pop()

    # 4. Start the search from every cell containing the start_letter
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == start_letter:
                # The initial call uses the trie_root, the function will handle the first letter
                find_words_recursive(r, c, trie_root, [])

    # 5. Determine and print the longest word found
    if not found_words:
        print(f"No words starting with '{start_letter}' were found.")
    else:
        longest_word = max(found_words, key=len)
        print(longest_word)

if __name__ == "__main__":
    solve_boggle()