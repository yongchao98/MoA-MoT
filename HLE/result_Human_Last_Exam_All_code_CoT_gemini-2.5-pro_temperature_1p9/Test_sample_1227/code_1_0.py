import numpy as np
import sys

def solve_puzzle():
    """
    Calculates the number of possible 8x8 crossword grids that satisfy a given set of rules
    by using a recursive backtracking algorithm.
    """

    # It's good practice to increase recursion limit for deep recursion,
    # though our depth of 18 is very shallow and doesn't require this.
    sys.setrecursionlimit(2000)

    class CrosswordGridCounter:
        """
        A class to encapsulate the grid and counting logic.
        """

        def __init__(self, size):
            if size % 2 != 0 or size < 4:
                raise ValueError("Size must be an even integer of at least 4.")
            self.size = size
            # Initialize the grid with all black squares (1). The border is now set.
            self.grid = np.ones((size, size), dtype=np.int8)
            self.count = 0
            
            # For 180-degree rotational symmetry on an N-sized grid with a black border,
            # we only need to decide the color for the cells in the top half of the inner (N-2)x(N-2) grid.
            inner_size = size - 2
            self.independent_cells = []
            for r_inner in range(inner_size // 2):
                for c_inner in range(inner_size):
                    # Convert inner grid coordinates to full grid coordinates
                    self.independent_cells.append((r_inner + 1, c_inner + 1))
            
            # For an 8x8 grid, this results in 3 rows * 6 cols = 18 independent cells to determine.

        def _is_connected(self):
            """Checks if all white squares (0s) form a single connected component using BFS."""
            white_squares = list(zip(*np.where(self.grid == 0)))
            if not white_squares:
                return False

            total_white_squares = len(white_squares)
            q = [white_squares[0]]
            visited = {white_squares[0]}
            
            head = 0
            while head < len(q):
                r, c = q[head]
                head += 1
                for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                    nr, nc = r + dr, c + dc
                    if self.grid[nr, nc] == 0 and (nr, nc) not in visited:
                        visited.add((nr, nc))
                        q.append((nr, nc))
            
            return len(visited) == total_white_squares

        def _has_valid_word_lengths(self):
            """Checks that all white square runs are at least 3 squares long."""
            for i in range(self.size):
                # Check row i
                length = 0
                for j in range(self.size):
                    if self.grid[i, j] == 0:
                        length += 1
                    else:
                        if 0 < length < 3: return False
                        length = 0
                if 0 < length < 3: return False

                # Check column i
                length = 0
                for j in range(self.size):
                    if self.grid[j, i] == 0:
                        length += 1
                    else:
                        if 0 < length < 3: return False
                        length = 0
                if 0 < length < 3: return False
            return True

        def _check_local_2x2_black_squares(self, r, c):
            """Checks if placing a black square at (r, c) creates a 2x2 black block."""
            for dr_start in [-1, 0]:
                for dc_start in [-1, 0]:
                    r_start, c_start = r + dr_start, c + dc_start
                    if (self.grid[r_start, c_start] == 1 and
                        self.grid[r_start + 1, c_start] == 1 and
                        self.grid[r_start, c_start + 1] == 1 and
                        self.grid[r_start + 1, c_start + 1] == 1):
                        return False
            return True

        def _solve(self, cell_index):
            """The recursive backtracking function."""
            if cell_index == len(self.independent_cells):
                # Base Case: A complete, symmetric grid is formed.
                # The 'no 2x2 black squares' check was already used for pruning.
                # Now perform the final connectivity and word length checks.
                if self._has_valid_word_lengths() and self._is_connected():
                    self.count += 1
                return

            r, c = self.independent_cells[cell_index]
            r_sym, c_sym = self.size - 1 - r, self.size - 1 - c

            # --- Choice 1: Place a WHITE square ---
            self.grid[r, c], self.grid[r_sym, c_sym] = 0, 0
            self._solve(cell_index + 1)

            # --- Choice 2: Place a BLACK square ---
            self.grid[r, c], self.grid[r_sym, c_sym] = 1, 1
            if self._check_local_2x2_black_squares(r, c) and \
               self._check_local_2x2_black_squares(r_sym, c_sym):
                self._solve(cell_index + 1)

        def find_grid_count(self):
            """Initializes and starts the search, returning the final count."""
            self._solve(0)
            return self.count

    grid_size = 8
    counter = CrosswordGridCounter(grid_size)
    number_of_grids = counter.find_grid_count()

    print(number_of_grids)

solve_puzzle()