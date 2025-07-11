import numpy as np
import sys

# Increase recursion limit for deep search, though good pruning should prevent hitting it.
sys.setrecursionlimit(2000)

class CrosswordGridEnumerator:
    """
    This class is designed to solve a specific combinatorial problem:
    counting the number of valid crossword puzzle grids under a set of rules.
    The rules for the 8x8 grid are:
    1. 180-degree rotational symmetry.
    2. A minimum word length of 3 (no runs of 1 or 2 white squares).
    3. Full interconnectivity of all white squares.
    4. No "cheater" squares, interpreted as no 2x2 blocks of black squares.
    
    The method used is a recursive backtracking search with aggressive pruning.
    """
    def __init__(self, size):
        if size % 2 != 0:
            raise ValueError("Size must be even for this implementation of symmetry.")
        self.size = size
        self.grid = np.full((size, size), -1, dtype=np.int8)
        self.count = 0
        
        # We only need to decide the state of the top half of the grid.
        # The bottom half is determined by symmetry.
        self.cells_to_set = []
        for r in range(self.size // 2):
            for c in range(self.size):
                self.cells_to_set.append((r, c))

    def solve(self):
        """
        Starts the backtracking search and returns the total count of valid grids.
        """
        self.backtrack(0)
        return self.count

    def backtrack(self, k):
        """
        The core recursive function for the backtracking search.
        'k' is the index of the independent cell we are currently deciding.
        """
        # Base case: If all independent cells are set, the grid is complete.
        if k == len(self.cells_to_set):
            # A full grid has been formed. Perform final validation.
            if self._is_fully_valid():
                self.count += 1
            return

        r, c = self.cells_to_set[k]
        r_sym, c_sym = self.size - 1 - r, self.size - 1 - c

        # --- Branch 1: Try placing a BLACK square (0) ---
        self.grid[r, c] = 0
        self.grid[r_sym, c_sym] = 0
        
        # Prune this branch if the new black square creates an immediate violation.
        if self._is_locally_valid(r, c):
            self.backtrack(k + 1)

        # --- Branch 2: Try placing a WHITE square (1) ---
        self.grid[r, c] = 1
        self.grid[r_sym, c_sym] = 1
        self.backtrack(k + 1)
        
    def _is_locally_valid(self, r, c):
        """
        Performs local checks to prune the search space after placing a BLACK square.
        """
        # Check for illegal short horizontal words closed by grid[r,c].
        if c > 0 and self.grid[r, c - 1] == 1:
            if c == 1 or self.grid[r, c - 2] == 0:  # ...B-W-B... or Edge-W-B
                return False
            if c > 1 and self.grid[r, c - 2] == 1 and (c == 2 or self.grid[r, c - 3] == 0): # ...B-WW-B... or Edge-WW-B
                return False

        # Check for illegal short vertical words.
        if r > 0 and self.grid[r - 1, c] == 1:
            if r == 1 or self.grid[r - 2, c] == 0:
                return False
            if r > 1 and self.grid[r - 2, c] == 1 and (r == 2 or self.grid[r - 3, c] == 0):
                return False

        # Check for 2x2 black blocks where (r,c) is the bottom-right corner.
        if r > 0 and c > 0 and self.grid[r - 1, c] == 0 and self.grid[r, c - 1] == 0 and self.grid[r - 1, c - 1] == 0:
            return False
        
        return True

    def _is_fully_valid(self):
        """
        Performs all validity checks on a completed grid.
        """
        if not self._is_structurally_sound():
            return False

        white_squares = list(zip(*np.where(self.grid == 1)))
        
        if not white_squares:
            return False

        q = [white_squares[0]]
        visited = {white_squares[0]}
        
        head = 0
        while head < len(q):
            curr_r, curr_c = q[head]
            head += 1
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = curr_r + dr, curr_c + dc
                if 0 <= nr < self.size and 0 <= nc < self.size and \
                   self.grid[nr, nc] == 1 and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))
                    
        return len(visited) == len(white_squares)

    def _is_structurally_sound(self):
        """
        Checks the entire grid for valid word lengths and no 2x2 black blocks.
        """
        for i in range(self.size):
            if any(0 < length < 3 for length in self._get_runs(self.grid[i, :])): return False
            if any(0 < length < 3 for length in self._get_runs(self.grid[:, i])): return False
        
        for r in range(self.size - 1):
            for c in range(self.size - 1):
                if np.all(self.grid[r:r+2, c:c+2] == 0):
                    return False
        return True

    def _get_runs(self, line):
        """Helper to find lengths of consecutive white squares (1s) in a line."""
        return [len(run) for run in ''.join(map(str, line)).split('0') if run]

if __name__ == '__main__':
    grid_size = 8
    print(f"Starting enumeration for {grid_size}x{grid_size} crossword grids...")
    print("This computation may take a few minutes...")
    
    enumerator = CrosswordGridEnumerator(grid_size)
    total_grids = enumerator.solve()

    print("\n--- Calculation Complete ---")
    print(f"The number of possible {grid_size}x{grid_size} grids is: {total_grids}")
<<<98>>>