import numpy as np

class CrosswordGridCounter:
    """
    A class to count the number of valid 8x8 crossword grids with specific constraints.
    """
    def __init__(self, size):
        """
        Initializes the grid counter for a given size.
        """
        if size % 2 != 0:
            raise ValueError("Grid size must be even for this symmetry implementation.")
        self.N = size
        self.grid = np.zeros((self.N, self.N), dtype=np.int8)
        self.count = 0

        # Assumption: The puzzle has a solid black border.
        self.grid[0, :] = 1
        self.grid[self.N - 1, :] = 1
        self.grid[:, 0] = 1
        self.grid[:, self.N - 1] = 1

        # Determine the cells we need to make decisions for.
        # Due to symmetry, this is the top half of the inner (N-2)x(N-2) grid.
        self.decision_cells = []
        inner_N = self.N - 2
        for r_inner in range(inner_N // 2):
            for c_inner in range(inner_N):
                # Map from inner grid coordinates to absolute grid coordinates
                self.decision_cells.append((r_inner + 1, c_inner + 1))
    
    def solve(self):
        """
        Starts the backtracking search and returns the total count of valid grids.
        """
        self.backtrack(0)
        return self.count

    def backtrack(self, k):
        """
        Recursively explores all possible grid configurations.
        """
        if k == len(self.decision_cells):
            # Base case: The grid is fully populated. Now, validate it.
            if self._is_valid():
                self.count += 1
            return

        r, c = self.decision_cells[k]
        r_sym, c_sym = self.N - 1 - r, self.N - 1 - c

        # Branch 1: Try placing a black square
        self.grid[r, c] = 1
        self.grid[r_sym, c_sym] = 1
        self.backtrack(k + 1)
        
        # Branch 2: Try placing a white square
        self.grid[r, c] = 0
        self.grid[r_sym, c_sym] = 0
        self.backtrack(k + 1)

    def _is_valid(self):
        """
        Checks if the current grid configuration is valid by applying all rules.
        """
        if not self._check_word_length():
            return False
        if not self._check_cheaters():
            return False
        # Connectivity is often the most expensive check, so do it last.
        if not self._check_connectivity():
            return False
        return True
    
    def _check_word_length(self):
        """
        Ensures no words (runs of white squares) are of length 1 or 2.
        """
        for grid_view in [self.grid, self.grid.T]: # Check rows, then columns via transpose
            for i in range(self.N):
                row_str = "".join(map(str, grid_view[i]))
                white_segments = row_str.split('1')
                if any(0 < len(segment) < 3 for segment in white_segments):
                    return False
        return True

    def _check_connectivity(self):
        """
        Ensures all white squares are connected into a single component.
        """
        white_squares = np.argwhere(self.grid == 0)
        if len(white_squares) == 0:
            # A grid with no words is invalid.
            return False

        start_node = tuple(white_squares[0])
        q = [start_node]
        visited = {start_node}
        head = 0
        while head < len(q):
            r, c = q[head]
            head += 1
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if self.grid[nr, nc] == 0 and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))
        
        return len(visited) == len(white_squares)

    def _check_cheaters(self):
        """
        Ensures there are no "cheater" black squares.
        This includes checking for 2x2 black blocks and isolated black squares.
        """
        for r in range(self.N - 1):
            for c in range(self.N - 1):
                # Rule 1: No 2x2 blocks of black squares.
                if self.grid[r, c] == 1 and self.grid[r+1, c] == 1 and self.grid[r, c+1] == 1 and self.grid[r+1, c+1] == 1:
                    return False
        
        for r in range(1, self.N - 1):
            for c in range(1, self.N - 1):
                if self.grid[r, c] == 1:
                    # Rule 2: No vertical cheaters (W on top/bottom, B on left/right)
                    if self.grid[r-1, c] == 0 and self.grid[r+1, c] == 0 and self.grid[r, c-1] == 1 and self.grid[r, c+1] == 1:
                        return False
                    # Rule 3: No horizontal cheaters (W on left/right, B on top/bottom)
                    if self.grid[r, c-1] == 0 and self.grid[r, c+1] == 0 and self.grid[r-1, c] == 1 and self.grid[r+1, c] == 1:
                        return False
        return True

if __name__ == '__main__':
    grid_finder = CrosswordGridCounter(8)
    total_grids = grid_finder.solve()
    print(total_grids)