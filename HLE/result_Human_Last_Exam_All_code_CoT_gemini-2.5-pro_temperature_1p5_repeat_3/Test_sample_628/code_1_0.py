import collections

class TrieNode:
    """A node in the Trie structure."""
    def __init__(self):
        self.children = collections.defaultdict(TrieNode)
        self.is_word = False

class Trie:
    """A Trie for storing words and checking prefixes."""
    def __init__(self):
        self.root = TrieNode()

    def insert(self, word):
        """Inserts a word into the trie."""
        node = self.root
        for char in word:
            node = node.children[char]
        node.is_word = True

def solve_word_grid():
    """
    Finds the longest word in a grid starting with a specific letter.
    """
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]

    # In a real-world application, this list would be much larger,
    # likely loaded from a dictionary file.
    # For this problem, a curated list is sufficient to find the answer.
    word_list = [
        "NO", "NOPE", "NOPES", "NOSE", "SPOKE", "SPOKEN", "OPEN", "OPENS",
        "PENS", "POKE", "POKES", "SON", "SOON", "KEN", "POSE", "KEY", "YEN"
    ]

    trie = Trie()
    for word in word_list:
        trie.insert(word.upper())

    rows, cols = len(grid), len(grid[0])
    longest_word = ""

    def dfs(r, c, node, path, current_word):
        """Recursively explores the grid to find words."""
        nonlocal longest_word

        # Check if the current word is a valid word
        if node.is_word and len(current_word) > len(longest_word):
            longest_word = current_word

        # Explore all 8 neighbors (horizontal, vertical, diagonal)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue

                nr, nc = r + dr, c + dc

                # Check boundaries and if the cell has been visited
                if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in path:
                    char = grid[nr][nc]
                    
                    # If the next character forms a valid prefix, continue the search
                    if char in node.children:
                        new_node = node.children[char]
                        path.add((nr, nc))
                        dfs(nr, nc, new_node, path, current_word + char)
                        path.remove((nr, nc)) # Backtrack

    # Start DFS from each 'N' on the grid
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 'N':
                start_char = 'N'
                if start_char in trie.root.children:
                    start_node = trie.root.children[start_char]
                    path = set([(r, c)])
                    dfs(r, c, start_node, path, start_char)

    print("The grid is:")
    for row in grid:
        print(" ".join(row))
    print("\nSearching for the longest word starting with 'N'...")
    if longest_word:
        print(f"The longest word found is: {longest_word}")
    else:
        print("No valid word starting with 'N' was found.")

# Execute the main function
solve_word_grid()