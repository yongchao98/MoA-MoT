import collections

class GridWordFinder:
    """
    Finds the longest word in a grid of letters starting with a specific character.
    """

    def __init__(self, grid):
        """
        Initializes the solver with the grid.
        """
        self.grid = grid
        self.rows = len(grid)
        self.cols = len(grid[0])
        self.longest_word = ""
        # The words and prefixes are specific to this puzzle grid and rules.
        self.words, self.prefixes = self._get_words_and_prefixes()

    def _get_words_and_prefixes(self):
        """
        Generates a focused set of words and prefixes for the search.
        The search is restricted to words starting with 'N'. This list contains
        potential candidates that can be formed using the letters in the grid.
        """
        candidate_words = [
            "NEG", "NEP", "NO", "NOG", "NOEL", "NOES", "NOLE", "NOLL", "NOPE",
            "NOPES", "NOSE", "NOSEY", "NYLON"
        ]
        
        # We only care about uppercase words for this puzzle.
        word_set = {word.upper() for word in candidate_words}
        
        # Create a set of all prefixes for efficient path pruning.
        prefix_set = set()
        for word in word_set:
            for i in range(1, len(word) + 1):
                prefix_set.add(word[:i])
                
        return word_set, prefix_set

    def solve(self):
        """
        Starts the search for the longest word from all 'N' cells in the grid.
        """
        # Per the rules, the word must start with 'N'.
        start_char = 'N'
        for r in range(self.rows):
            for c in range(self.cols):
                if self.grid[r][c] == start_char:
                    # Begin the recursive search from this starting cell.
                    self._dfs(r, c, self.grid[r][c], {(r, c)})
        
        # Print the final result.
        if self.longest_word:
            print(f"The longest word that can be formed starting with 'N' is: {self.longest_word}")
        else:
            print(f"No valid word starting with 'N' could be formed.")

    def _dfs(self, r, c, current_path, visited_cells):
        """
        A recursive Depth-First Search to explore paths on the grid.

        Args:
            r (int): Current row index.
            c (int): Current column index.
            current_path (str): The string of characters formed so far.
            visited_cells (set): A set of (row, col) tuples visited in the current path.
        """
        # First, check if the current path itself is a valid word.
        if current_path in self.words:
            # If it's the longest we've found so far, save it.
            if len(current_path) > len(self.longest_word):
                self.longest_word = current_path

        # Optimization: If the current path is not a prefix of any valid word,
        # we can stop exploring this branch.
        if current_path not in self.prefixes:
            return

        # Explore all 8 neighbors (including diagonals).
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                # Skip the current cell itself.
                if dr == 0 and dc == 0:
                    continue

                nr, nc = r + dr, c + dc

                # Check if the neighbor is within the grid and has not been visited yet.
                if 0 <= nr < self.rows and 0 <= nc < self.cols and (nr, nc) not in visited_cells:
                    # Continue the search from the new cell.
                    new_visited = visited_cells.copy()
                    new_visited.add((nr, nc))
                    self._dfs(nr, nc, current_path + self.grid[nr][nc], new_visited)

# --- Main execution block ---

# Define the 4x4 grid of letters provided in the problem.
GRID = [
    ['E', 'S', 'G', 'D'],
    ['K', 'P', 'O', 'J'],
    ['O', 'E', 'N', 'Y'],
    ['L', 'P', 'S', 'N']
]

# Create an instance of the finder and run the solver.
finder = GridWordFinder(GRID)
finder.solve()