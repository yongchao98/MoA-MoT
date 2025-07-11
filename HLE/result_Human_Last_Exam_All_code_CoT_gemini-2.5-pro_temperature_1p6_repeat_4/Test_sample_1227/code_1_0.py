import sys

class CrosswordGridCounter:
    """
    Solves the problem of counting valid 8x8 crossword grids with specific rules.

    This class implements a backtracking search assuming a standard crossword convention
    where the 8x8 grid has a black perimeter, reducing the problem to the inner 6x6 grid.
    This makes the search space computationally feasible.
    """

    def __init__(self):
        # Grid dimensions
        self.size = 8
        self.inner_size = 6
        
        # Initialize the 8x8 grid with a black perimeter
        self.grid = [[0] * self.size for _ in range(self.size)]
        
        # Counter for valid grids found
        self.valid_grid_count = 0

    def solve(self):
        """
        Starts the process of finding and counting all valid grid patterns.
        The result is printed to the console.
        """
        # The recursive backtracking starts from the top-left of the inner 6x6 grid.
        # The coordinates (r_inner, c_inner) are 0-indexed relative to the 6x6 grid.
        self._backtrack(0, 0)
        print(self.valid_grid_count)

    def _backtrack(self, r_inner, c_inner):
        """
        Recursively explores all symmetric patterns for the inner grid.
        """
        # Base case: We have filled the top half of the inner grid (3 rows).
        if r_inner >= self.inner_size // 2:
            self._validate_and_count()
            return

        # Map inner grid coordinates (0-5) to the full 8x8 grid coordinates (1-6)
        r_abs, c_abs = r_inner + 1, c_inner + 1

        # Determine the coordinates for the next recursive step
        next_r_inner, next_c_inner = (r_inner, c_inner + 1) if c_inner + 1 < self.inner_size else (r_inner + 1, 0)

        # --- Choice 1: Place a black square (0) ---
        self.grid[r_abs][c_abs] = 0
        self.grid[self.size - 1 - r_abs][self.size - 1 - c_abs] = 0
        self._backtrack(next_r_inner, next_c_inner)

        # --- Choice 2: Place a white square (1) ---
        self.grid[r_abs][c_abs] = 1
        self.grid[self.size - 1 - r_abs][self.size - 1 - c_abs] = 1
        self._backtrack(next_r_inner, next_c_inner)

    def _validate_and_count(self):
        """
        Performs all validation checks on the fully generated grid.
        """
        if (self._check_word_lengths() and
            self._check_connectivity() and
            not self._has_cheater_squares()):
            self.valid_grid_count += 1

    def _check_word_lengths(self):
        """Checks if all words are at least 3 letters long."""
        for r in range(self.size):
            for c in range(self.size):
                if self.grid[r][c] == 1:
                    # Check for the start of a horizontal word
                    if c == 0 or self.grid[r][c-1] == 0:
                        length = 0
                        while c + length < self.size and self.grid[r][c+length] == 1:
                            length += 1
                        if length < 3:
                            return False
                    # Check for the start of a vertical word
                    if r == 0 or self.grid[r-1][c] == 0:
                        length = 0
                        while r + length < self.size and self.grid[r+length][c] == 1:
                            length += 1
                        if length < 3:
                            return False
        return True

    def _check_connectivity(self):
        """Checks if all white squares are connected."""
        total_white = sum(row.count(1) for row in self.grid)
        if total_white == 0:
            return False  # A grid with no words is not a valid puzzle.

        q = []
        visited = set()
        
        # Find the first white square to start the search from
        for r_start in range(self.size):
            try:
                c_start = self.grid[r_start].index(1)
                q.append((r_start, c_start))
                visited.add((r_start, c_start))
                break
            except ValueError:
                continue
        
        # Perform BFS to find all reachable white squares
        head = 0
        while head < len(q):
            r, c = q[head]
            head += 1
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < self.size and 0 <= nc < self.size and \
                   self.grid[nr][nc] == 1 and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))
                    
        return len(visited) == total_white

    def _has_cheater_squares(self):
        """Checks for any "cheater" black squares."""
        # Cheater squares can only exist in the inner part of the grid
        for r in range(1, self.size - 1):
            for c in range(1, self.size - 1):
                if self.grid[r][c] == 0:  # It's a black square
                    # A square is NOT a cheater if it separates words.
                    separates_vertically = (self.grid[r-1][c] == 1 and self.grid[r+1][c] == 1)
                    separates_horizontally = (self.grid[r][c-1] == 1 and self.grid[r][c+1] == 1)
                    
                    if not separates_vertically and not separates_horizontally:
                        return True  # Found a cheater square
        return False

if __name__ == '__main__':
    counter = CrosswordGridCounter()
    counter.solve()