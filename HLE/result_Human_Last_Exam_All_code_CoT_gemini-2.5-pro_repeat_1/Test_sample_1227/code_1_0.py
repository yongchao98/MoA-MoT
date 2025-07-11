import sys

class CrosswordGridCounter:
    """
    This class calculates the number of possible 8x8 crossword grids
    that satisfy a specific set of rules:
    1. 180-degree rotational symmetry.
    2. Minimum word length of 3 (for both horizontal and vertical words).
    3. All white squares must be fully interconnected.
    4. No "cheater" squares, interpreted as no 2x2 blocks of black squares.
    """

    def __init__(self, size=8):
        """
        Initializes the grid counter.
        """
        if size % 2 != 0:
            raise ValueError("Grid size must be even for this symmetry implementation.")
        self.size = size
        self.grid = [[0] * size for _ in range(size)]  # 0 for white, 1 for black
        self.count = 0
        # For an N*N grid with 180-degree symmetry, we only need to decide
        # the state of the first N*N/2 cells. The rest are determined.
        self.independent_cells = (size * size) // 2

    def solve(self):
        """
        Starts the backtracking search and returns the total count of valid grids.
        """
        self.backtrack(0)
        return self.count

    def backtrack(self, index):
        """
        Recursively fills the grid, checking for violations at each step to prune
        the search space.
        
        Args:
            index (int): The index of the independent cell to be filled (0 to 31).
        """
        # Base case: All independent cells have been filled.
        if index == self.independent_cells:
            # The top half of the grid is filled, and the bottom half is determined by symmetry.
            # At this point, 2x2 and row-length checks have already passed during pruning.
            # We only need to perform the final checks on columns and connectivity.
            if self.check_column_word_lengths() and self.check_connectivity():
                self.count += 1
            return

        r = index // self.size
        c = index % self.size
        sym_r, sym_c = self.size - 1 - r, self.size - 1 - c

        # --- Choice 1: Place a BLACK square ---
        self.grid[r][c] = 1
        self.grid[sym_r][sym_c] = 1

        # Pruning check: See if this new black square creates a 2x2 block.
        # We only need to check the 2x2 square for which (r, c) is the bottom-right corner,
        # as the other cells in that potential square would have been filled previously.
        if not (r > 0 and c > 0 and self.grid[r - 1][c - 1] == 1 and
                self.grid[r - 1][c] == 1 and self.grid[r][c - 1] == 1):
            
            # Pruning check: If a row in the top half is now complete, check its word lengths.
            if c == self.size - 1:
                if self.check_single_row_length(r) and self.check_single_row_length(sym_r):
                    self.backtrack(index + 1)
            else:
                self.backtrack(index + 1)

        # --- Choice 2: Place a WHITE square (and backtrack from the black choice) ---
        self.grid[r][c] = 0
        self.grid[sym_r][sym_c] = 0

        # Pruning check: If a row is now complete, check its word lengths.
        if c == self.size - 1:
            if self.check_single_row_length(r) and self.check_single_row_length(sym_r):
                self.backtrack(index + 1)
        else:
            self.backtrack(index + 1)
        
    def check_single_row_length(self, r):
        """Checks that all horizontal words in a given row are of length 3 or more."""
        length = 0
        for c in range(self.size):
            if self.grid[r][c] == 0:  # White square
                length += 1
            else:  # Black square
                if 0 < length < 3:
                    return False
                length = 0
        return not (0 < length < 3)

    def check_column_word_lengths(self):
        """Checks that all vertical words in the entire grid are of length 3 or more."""
        for c in range(self.size):
            length = 0
            for r in range(self.size):
                if self.grid[r][c] == 0:  # White square
                    length += 1
                else:  # Black square
                    if 0 < length < 3:
                        return False
                    length = 0
            if 0 < length < 3:
                return False
        return True

    def check_connectivity(self):
        """Checks if all white squares form a single connected component using BFS."""
        white_squares = []
        first_white = None
        for r in range(self.size):
            for c in range(self.size):
                if self.grid[r][c] == 0:
                    if first_white is None:
                        first_white = (r, c)
                    white_squares.append((r, c))

        # If there are no white squares, it's not a valid puzzle grid.
        if not first_white:
            return False

        q = [first_white]
        visited = {first_white}
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

if __name__ == "__main__":
    counter = CrosswordGridCounter(size=8)
    total_grids = counter.solve()
    print(total_grids)