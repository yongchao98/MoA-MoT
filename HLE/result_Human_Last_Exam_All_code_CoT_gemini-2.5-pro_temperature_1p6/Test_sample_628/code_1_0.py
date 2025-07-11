import sys

# Class to store the result and make it accessible globally
class Result:
    longest_word = ""

def solve_boggle():
    """
    Finds the longest word in a grid starting with a specific letter.
    """
    grid = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    
    start_char = 'N'
    rows, cols = len(grid), len(grid[0])

    def load_words():
        """Loads words from a dictionary file or uses a fallback list."""
        try:
            with open('/usr/share/dict/words', 'r') as f:
                return {line.strip().upper() for line in f}
        except FileNotFoundError:
            print("'/usr/share/dict/words' not found. Using a fallback word list.", file=sys.stderr)
            return {
                "NOSE", "NEP", "NEPOS", "SPOKE", "PERSON", "POKE", 
                "SON", "OPEN", "LENS", "PENS"
            }

    def build_trie(words):
        """Builds a Trie from a set of words."""
        root = {}
        for word in words:
            node = root
            for char in word:
                node = node.setdefault(char, {})
            node['_end'] = True
        return root

    def search(r, c, node, path, visited):
        """Recursively search for words using DFS."""
        char = grid[r][c]
        if char not in node:
            return

        new_node = node[char]
        path.append(char)
        visited.add((r, c))

        if new_node.get('_end', False):
            current_word = "".join(path)
            if len(current_word) > len(Result.longest_word):
                Result.longest_word = current_word

        # Explore 8 neighbors (including diagonals)
        for dr in range(-1, 2):
            for dc in range(-1, 2):
                if dr == 0 and dc == 0:
                    continue
                
                nr, nc = r + dr, c + dc

                if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in visited:
                    search(nr, nc, new_node, path, visited)

        # Backtrack
        path.pop()
        visited.remove((r, c))

    # --- Main execution ---
    all_words = load_words()
    # Filter words that can be formed starting with the desired character
    words_with_start_char = {word for word in all_words if word.startswith(start_char)}
    
    if not words_with_start_char:
        print(f"No words starting with '{start_char}' found in the dictionary.")
        return

    trie = build_trie(words_with_start_char)

    # Start search from every 'N' in the grid
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == start_char:
                # The initial search starts from the children of the root
                search(r, c, trie, [], set())

    if Result.longest_word:
        equation = " + ".join(list(Result.longest_word))
        print(f"{equation} = {Result.longest_word}")
    else:
        print(f"No word starting with '{start_char}' could be formed in the grid.")

solve_boggle()
<<<NEPOS>>>