import sys

class WordFinder:
    """
    Finds the longest word in a grid of letters starting with a specific character.
    """

    def __init__(self, grid, start_char):
        self.grid = [[char.lower() for char in row] for row in grid]
        self.start_char = start_char.lower()
        self.rows = len(grid)
        self.cols = len(grid[0])
        self.found_words = set()
        self.word_set, self.prefix_set = self._load_dictionary()

    def _load_dictionary(self):
        """
        Loads words and prefixes from a dictionary file.
        Uses a fallback list if the file is not found.
        """
        word_set = set()
        prefix_set = set()
        try:
            with open('/usr/share/dict/words', 'r') as f:
                for word in f:
                    # We only care about words that could be formed from the grid
                    w = word.strip().lower()
                    if len(w) >= 3:  # At least 3 letters for a meaningful word
                        word_set.add(w)
                        for i in range(1, len(w) + 1):
                            prefix_set.add(w[:i])
        except FileNotFoundError:
            sys.stderr.write("Warning: Dictionary file '/usr/share/dict/words' not found.\n")
            sys.stderr.write("Using a smaller, built-in word list.\n")
            fallback_words = ['esg', 'esk', 'esp', 'ego', 'eon', 'son', 'sony', 'song',
                              'pes', 'peg', 'pen', 'peon', 'peso', 'pos', 'pose', 'pons',
                              'goes', 'gone', 'gen', 'yen', 'yep', 'yens', 'keps',
                              'nope', 'nopes', 'nose', 'nosy', 'sonly', 'syne']
            for w in fallback_words:
                word_set.add(w)
                for i in range(1, len(w) + 1):
                    prefix_set.add(w[:i])

        # Add all prefixes for words starting with the required letter
        final_prefix_set = {p for p in prefix_set if p.startswith(self.start_char)}
        final_word_set = {w for w in word_set if w.startswith(self.start_char)}
        return final_word_set, final_prefix_set
    
    def _is_valid(self, r, c, visited):
        """Checks if a cell is within grid bounds and hasn't been visited."""
        return 0 <= r < self.rows and 0 <= c < self.cols and (r, c) not in visited

    def _search(self, r, c, current_word, visited):
        """
        Performs the recursive Depth-First Search.
        """
        # Append the new character to the current word
        current_word += self.grid[r][c]

        # Prune the search if the current path is not a prefix of any word
        if current_word not in self.prefix_set:
            return

        # If the current path forms a valid word, add it to our findings
        if current_word in self.word_set:
            self.found_words.add(current_word)
        
        # Mark the current cell as visited
        visited.add((r, c))

        # Explore all 8 neighbors (horizontal, vertical, diagonal)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue  # Skip the cell itself
                
                next_r, next_c = r + dr, c + dc

                if self._is_valid(next_r, next_c, visited):
                    self._search(next_r, next_c, current_word, visited.copy())

    def solve(self):
        """
        Initiates the search from all valid starting points and prints the result.
        """
        # Find all starting positions for the given character
        start_positions = []
        for r in range(self.rows):
            for c in range(self.cols):
                if self.grid[r][c] == self.start_char:
                    start_positions.append((r, c))
        
        if not start_positions:
            print(f"Starting letter '{self.start_char.upper()}' not found in grid.")
            return
            
        # Start the DFS from each starting position
        for r_start, c_start in start_positions:
            self._search(r_start, c_start, "", set())
            
        # Find the longest word among all that were found
        if not self.found_words:
            print("No words starting with 'N' found.")
        else:
            longest_word = ""
            for word in self.found_words:
                if len(word) > len(longest_word):
                    longest_word = word
            
            # Print each letter of the final word as per instruction format.
            print("The longest word found is:", longest_word)
            print("The letters are:", " ".join(list(longest_word.upper())))


if __name__ == "__main__":
    grid_data = [
        ['E', 'S', 'G', 'D'],
        ['K', 'P', 'O', 'J'],
        ['O', 'E', 'N', 'Y'],
        ['L', 'P', 'S', 'N']
    ]
    
    finder = WordFinder(grid=grid_data, start_char='N')
    finder.solve()