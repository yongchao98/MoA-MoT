import sys

# It is a good practice to increase the recursion limit for deep search problems.
sys.setrecursionlimit(2000)

class CrosswordGridCounter:
    """
    This class finds the number of valid crossword grids for a given size and a specific set of rules
    by using a recursive backtracking algorithm.
    """

    def __init__(self, size):
        """Initializes the solver for a grid of a given size."""
        if size % 2 != 0:
            raise ValueError("Grid size must be even for this implementation of 180-degree symmetry.")
        self.N = size
        self.grid = [[-1] * self.N for _ in range(self.N)] # -1 indicates an undecided square
        self.solution_count = 0

        # We only need to decide the color for cells in one half of the grid due to symmetry.
        # For an 8x8 grid with pre-filled borders, this leaves 18 cells to determine.
        self.cells_to_fill = []
        for r in range(1, self.N // 2):
            for c in range(1, self.N - 1):
                self.cells_to_fill.append((r, c))

    def solve(self):
        """Starts the recursive search and returns the total count of valid grids."""
        # A standard crossword grid convention is to have the border be all black squares.
        # We use 1 for black and 0 for white.
        for i in range(self.N):
            self.grid[0][i] = 1
            self.grid[self.N - 1][i] = 1
            self.grid[i][0] = 1
            self.grid[i][self.N - 1] = 1
        
        self.generate(0)
        return self.solution_count

    def generate(self, cell_index):
        """
        Recursively explores possibilities by placing white or black squares.
        This function constitutes the core of the backtracking search.
        """
        # Base case: If all independent cells have been filled, the grid is complete.
        # Now we perform the final checks.
        if cell_index == len(self.cells_to_fill):
            if self._is_grid_valid():
                self.solution_count += 1
            return

        r, c = self.cells_to_fill[cell_index]
        r_sym = self.N - 1 - r
        c_sym = self.N - 1 - c

        # --- Branch 1: Try placing a white square ---
        self.grid[r][c] = 0
        self.grid[r_sym][c_sym] = 0
        self.generate(cell_index + 1)

        # --- Branch 2: Try placing a black square ---
        self.grid[r][c] = 1
        self.grid[r_sym][c_sym] = 1
        # Pruning Step: We check for the "no 2x2 black squares" rule immediately.
        # If this rule is violated, we don't need to explore this path any further.
        if not self._has_2x2_violation(r, c) and not self._has_2x2_violation(r_sym, c_sym):
            self.generate(cell_index + 1)

    def _is_grid_valid(self):
        """Checks if a fully generated grid meets all final criteria."""
        return self._check_word_lengths() and self._check_connectivity()

    def _check_word_lengths(self):
        """Ensures all words (runs of white squares) are at least 3 squares long."""
        # Check all rows for sequences of 1 or 2 white squares
        for r in range(self.N):
            length = 0
            for c in range(self.N):
                if self.grid[r][c] == 0:
                    length += 1
                else:
                    if 0 < length < 3: return False
                    length = 0
            if 0 < length < 3: return False
        
        # Check all columns for sequences of 1 or 2 white squares
        for c in range(self.N):
            length = 0
            for r in range(self.N):
                if self.grid[r][c] == 0:
                    length += 1
                else:
                    if 0 < length < 3: return False
                    length = 0
            if 0 < length < 3: return False
        return True

    def _check_connectivity(self):
        """Ensures all white squares form a single connected component using BFS."""
        first_white = None
        total_white = 0
        for r in range(self.N):
            for c in range(self.N):
                if self.grid[r][c] == 0:
                    if first_white is None:
                        first_white = (r, c)
                    total_white += 1
        
        # A valid puzzle must have empty spaces for words. An all-black grid is not a puzzle.
        if total_white == 0:
            return False

        # Use BFS to find all connected white squares from the first one found.
        q = [first_white]
        visited = {first_white}
        connected_count = 1
        head = 0
        while head < len(q):
            r, c = q[head]
            head += 1
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < self.N and 0 <= nc < self.N and self.grid[nr][nc] == 0 and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))
                    connected_count += 1
        
        # The grid is fully connected if the number of visited squares equals the total number.
        return connected_count == total_white

    def _has_2x2_violation(self, r_new, c_new):
        """Checks if placing a black square at (r_new, c_new) creates a 2x2 block."""
        # A new black square can complete a 2x2 block in four ways (by being any of its four corners).
        # Check if the new square is the bottom-right corner of a 2x2 block
        if r_new > 0 and c_new > 0 and \
           self.grid[r_new - 1][c_new] == 1 and self.grid[r_new][c_new - 1] == 1 and self.grid[r_new - 1][c_new - 1] == 1:
            return True
        # Check if it's the bottom-left corner
        if r_new > 0 and c_new < self.N - 1 and \
           self.grid[r_new - 1][c_new] == 1 and self.grid[r_new][c_new + 1] == 1 and self.grid[r_new - 1][c_new + 1] == 1:
            return True
        # Check if it's the top-right corner
        if r_new < self.N - 1 and c_new > 0 and \
           self.grid[r_new + 1][c_new] == 1 and self.grid[r_new][c_new - 1] == 1 and self.grid[r_new + 1][c_new - 1] == 1:
            return True
        # Check if it's the top-left corner
        if r_new < self.N - 1 and c_new < self.N - 1 and \
           self.grid[r_new + 1][c_new] == 1 and self.grid[r_new][c_new + 1] == 1 and self.grid[r_new + 1][c_new + 1] == 1:
            return True
        return False

# This block allows the script to be executed directly to find the solution.
if __name__ == '__main__':
    grid_size = 8
    solver = CrosswordGridCounter(grid_size)
    number_of_grids = solver.solve()
    print(f"The number of possible {grid_size}x{grid_size} crossword grids is: {number_of_grids}")
