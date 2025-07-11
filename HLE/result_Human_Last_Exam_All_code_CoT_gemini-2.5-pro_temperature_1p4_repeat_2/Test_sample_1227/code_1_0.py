import sys

# We increase the recursion limit for the backtracking algorithm.
# The default limit might be too low for the depth of this search.
sys.setrecursionlimit(2000)

class CrosswordGridCounter:
    """
    Solves the crossword grid counting problem using a pruned backtracking algorithm.
    """
    def __init__(self, size):
        if size < 3 or size % 2 != 0:
            # The described symmetry works best for even-sized grids.
            # Odd-sized grids that meet these rules are rare (e.g., none for 3, 5, 7).
            raise ValueError("Size must be an even integer >= 4.")
        
        self.N = size
        # Grid state: -1 for undecided, 0 for black, 1 for white.
        self.grid = [[-1] * self.N for _ in range(self.N)]
        self.count = 0
        # List of independent cells to decide due to 180-degree symmetry.
        self.cells_to_fill = [(r, c) for r in range(self.N // 2) for c in range(self.N)]

    def solve_and_print(self):
        """Starts the solver and prints the final result."""
        self._backtrack(0)
        print(f"The number of possible {self.N}x{self.N} grids is: {self.count}")

    def _backtrack(self, k):
        """
        Recursively explores grid configurations.
        'k' is the index of the cell in 'cells_to_fill' we are currently deciding.
        """
        # Base Case: All independent cells have been decided.
        if k == len(self.cells_to_fill):
            if self._is_valid_final_grid():
                self.count += 1
            return

        r, c = self.cells_to_fill[k]
        sr, sc = self.N - 1 - r, self.N - 1 - c

        # --- Branch 1: Try placing a black square ---
        self.grid[r][c] = self.grid[sr][sc] = 0
        # Prune if this placement creates a 2x2 block of black squares.
        if not self._has_cheater_violation(r, c) and not self._has_cheater_violation(sr, sc):
            # If a row has just been completed, prune if it has short words.
            if c == self.N - 1:
                if not self._has_row_short_word_violation(r) and not self._has_row_short_word_violation(sr):
                    self._backtrack(k + 1)
            else:
                self._backtrack(k + 1)

        # --- Branch 2: Try placing a white square ---
        self.grid[r][c] = self.grid[sr][sc] = 1
        # Prune if the row just completed has short words.
        if c == self.N - 1:
            if not self._has_row_short_word_violation(r) and not self._has_row_short_word_violation(sr):
                self._backtrack(k + 1)
        else:
            self._backtrack(k + 1)
        
        # Backtrack: Reset the cells to undecided for the parent recursion call.
        self.grid[r][c] = self.grid[sr][sc] = -1

    def _has_cheater_violation(self, r, c):
        """Checks if placing a black square at (r, c) creates a 2x2 black block."""
        for dr in [-1, 0]:
            for dc in [-1, 0]:
                r_start, c_start = r + dr, c + dc
                if 0 <= r_start < self.N - 1 and 0 <= c_start < self.N - 1:
                    if self.grid[r_start][c_start] == 0 and \
                       self.grid[r_start + 1][c_start] == 0 and \
                       self.grid[r_start][c_start + 1] == 0 and \
                       self.grid[r_start + 1][c_start + 1] == 0:
                        return True
        return False

    def _has_row_short_word_violation(self, r):
        """Checks a fully decided row for any white segments of length 1 or 2."""
        length = 0
        for c in range(self.N):
            if self.grid[r][c] == 1:
                length += 1
            else:
                if 0 < length < 3: return True
                length = 0
        return 0 < length < 3

    def _is_valid_final_grid(self):
        """Performs final validation for columns and connectivity on a complete grid."""
        # 1. Check columns for short words.
        for c in range(self.N):
            length = 0
            for r in range(self.N):
                if self.grid[r][c] == 1:
                    length += 1
                else:
                    if 0 < length < 3: return False
                    length = 0
            if 0 < length < 3: return False
        
        # 2. Check for connectivity of white squares.
        white_squares = [(r, c) for r in range(self.N) for c in range(self.N) if self.grid[r][c] == 1]
        if not white_squares:
            return False  # No words are not allowed.

        # BFS/flood-fill to check connectivity.
        q = [white_squares[0]]
        visited = {white_squares[0]}
        head = 0
        while head < len(q):
            r_node, c_node = q[head]
            head += 1
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r_node + dr, c_node + dc
                if 0 <= nr < self.N and 0 <= nc < self.N and self.grid[nr][nc] == 1 and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))
        
        return len(visited) == len(white_squares)

if __name__ == '__main__':
    # Initialize the counter for an 8x8 grid.
    counter = CrosswordGridCounter(size=8)
    # Run the solver.
    counter.solve_and_print()