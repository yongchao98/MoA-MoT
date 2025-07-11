import sys

class CrosswordGridCounter:
    """
    Solves for the number of valid crossword grids under a specific set of rules.
    """
    def __init__(self, size):
        """
        Initializes the solver for a grid of a given size.
        - size: The width and height of the grid (e.g., 8 for 8x8).
        - grid: Represents the puzzle, with 1 for black, 0 for white, -1 for undecided.
        - solution_count: The final count of valid grids.
        - cells_to_decide: The number of unique cells that determine the grid due to symmetry.
        """
        if size % 2 != 0:
            raise ValueError("Grid size must be even for this symmetry implementation.")
        self.size = size
        self.grid = [[-1] * size for _ in range(size)]
        self.solution_count = 0
        self.cells_to_decide = (size * size) // 2

    def count_grids(self):
        """
        Public method to start the recursive search and return the result.
        """
        self._recurse(0)
        return self.solution_count

    def _recurse(self, index):
        """
        Recursively places squares to build and validate all possible grids.
        'index' corresponds to a cell in the top half of the grid (0 to 31 for 8x8).
        """
        # Base case: All independent cells have been assigned a color.
        # Now, perform final validation on the complete grid.
        if index == self.cells_to_decide:
            if self._is_valid_final_grid():
                self.solution_count += 1
            return

        # Convert the linear index to 2D coordinates for the top half of the grid.
        r = index // self.size
        c = index % self.size
        
        # Determine the coordinates of the symmetrically opposite square.
        sr, sc = self.size - 1 - r, self.size - 1 - c

        # --- Choice 1: Place a black square ---
        self.grid[r][c] = 1
        self.grid[sr][sc] = 1

        # Pruning Step: If placing this black square creates a 2x2 block of black
        # squares (a "cheater" square), abandon this path.
        if not self._has_2x2_black(r, c) and not self._has_2x2_black(sr, sc):
            self._recurse(index + 1)

        # --- Choice 2: Place a white square ---
        self.grid[r][c] = 0
        self.grid[sr][sc] = 0
        self._recurse(index + 1)

    def _has_2x2_black(self, r, c):
        """
        Checks if the newly placed black square at (r, c) completes any 2x2 block.
        It checks the four 2x2 squares that (r, c) could be a corner of.
        """
        for dr in [-1, 0]:
            for dc in [-1, 0]:
                r_start, c_start = r + dr, c + dc
                if 0 <= r_start < self.size - 1 and 0 <= c_start < self.size - 1:
                    if (self.grid[r_start][c_start] == 1 and
                        self.grid[r_start + 1][c_start] == 1 and
                        self.grid[r_start][c_start + 1] == 1 and
                        self.grid[r_start + 1][c_start + 1] == 1):
                        return True
        return False

    def _is_valid_final_grid(self):
        """
        Performs the final, more expensive checks on a fully generated grid.
        """
        return self._check_word_lengths() and self._check_connectivity()

    def _check_word_lengths(self):
        """
        Checks if all words (runs of white squares) are of length 3 or more.
        It pads rows and columns with black squares to handle edges correctly.
        """
        for i in range(self.size):
            # Check row i
            row_vals = [self.grid[i][c] for c in range(self.size)]
            if self._has_short_words_in_line(row_vals):
                return False
            
            # Check column i
            col_vals = [self.grid[r][i] for r in range(self.size)]
            if self._has_short_words_in_line(col_vals):
                return False
        return True

    def _has_short_words_in_line(self, line):
        """Helper for _check_word_lengths. Checks a single row or column."""
        padded_line = [1] + line + [1]
        for i in range(len(padded_line) - 3):
            # Check for BWB (word of length 1) -> [1, 0, 1]
            if padded_line[i] == 1 and padded_line[i+1] == 0 and padded_line[i+2] == 1:
                return True
            # Check for BWWB (word of length 2) -> [1, 0, 0, 1]
            if padded_line[i] == 1 and padded_line[i+1] == 0 and padded_line[i+2] == 0 and padded_line[i+3] == 1:
                return True
        return False

    def _check_connectivity(self):
        """
        Checks if all white squares form a single contiguous block using BFS.
        """
        white_squares = []
        for r in range(self.size):
            for c in range(self.size):
                if self.grid[r][c] == 0:
                    white_squares.append((r, c))
        
        # If there are no white squares, it's not a valid puzzle grid.
        if not white_squares:
            return False

        q = [white_squares[0]]
        visited = {white_squares[0]}
        
        head = 0
        while head < len(q):
            r, c = q[head]
            head += 1
            
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if (0 <= nr < self.size and 0 <= nc < self.size and
                        self.grid[nr][nc] == 0 and (nr, nc) not in visited):
                    visited.add((nr, nc))
                    q.append((nr, nc))

        return len(visited) == len(white_squares)

if __name__ == '__main__':
    # Note: This is a computationally intensive task.
    # The search space is 2^32, and while pruning helps, this will take a long time to run.
    grid_size = 8
    counter = CrosswordGridCounter(grid_size)
    total_grids = counter.count_grids()
    print(total_grids)
