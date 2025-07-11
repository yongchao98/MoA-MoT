import sys

class CrosswordGridCounter:
    """
    This class calculates the number of valid 8x8 crossword grids
    based on a set of constraints, assuming a standard black border.
    """
    def __init__(self, size=8):
        self.size = size
        self.grid = [['' for _ in range(size)] for _ in range(size)]
        self.valid_grids_count = 0

        # Assumption: Standard American-style crosswords have a fully black border.
        # This reduces the problem space from 2^32 to a tractable 2^18.
        for i in range(self.size):
            self.grid[i][0] = 'B'
            self.grid[i][self.size - 1] = 'B'
            self.grid[0][i] = 'B'
            self.grid[self.size - 1][i] = 'B'

        # Identify the unique cells that determine the grid pattern due to symmetry.
        # For an 8x8 grid with black borders, we decide the top half of the inner 6x6 grid.
        # These are rows 1, 2, 3 and columns 1 through 6, totaling 18 cells.
        self.cells_to_fill = []
        for r in range(1, self.size // 2):
            for c in range(1, self.size - 1):
                self.cells_to_fill.append((r, c))
        
    def _is_valid_coord(self, r, c):
        return 0 <= r < self.size and 0 <= c < self.size

    def _get_color(self, r, c):
        """Helper to get a grid color, treating boundaries as Black."""
        if self._is_valid_coord(r, c):
            return self.grid[r][c]
        return 'B'

    def check_word_lengths(self):
        """Ensures all words (runs of white squares) are at least 3 letters long."""
        for r in range(self.size):
            length = 0
            for c in range(self.size):
                if self.grid[r][c] == 'W':
                    length += 1
                else:
                    if 1 <= length <= 2:
                        return False
                    length = 0
            if 1 <= length <= 2:
                return False
        
        for c in range(self.size):
            length = 0
            for r in range(self.size):
                if self.grid[r][c] == 'W':
                    length += 1
                else:
                    if 1 <= length <= 2:
                        return False
                    length = 0
            if 1 <= length <= 2:
                return False
        return True

    def check_connectivity(self):
        """Ensures all white squares form a single connected component."""
        white_squares = []
        for r in range(self.size):
            for c in range(self.size):
                if self.grid[r][c] == 'W':
                    white_squares.append((r, c))
        
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
                if self._is_valid_coord(nr, nc) and self.grid[nr][nc] == 'W' and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))
        
        return len(visited) == len(white_squares)
        
    def check_no_cheaters(self):
        """Ensures no black squares form thin walls between white squares."""
        for r in range(self.size):
            for c in range(self.size):
                if self.grid[r][c] == 'B':
                    is_vert_wall = (self._get_color(r - 1, c) == 'W' and self._get_color(r + 1, c) == 'W')
                    is_horiz_blocked = (self._get_color(r, c - 1) == 'B' and self._get_color(r, c + 1) == 'B')
                    if is_vert_wall and is_horiz_blocked:
                        return False
                        
                    is_horiz_wall = (self._get_color(r, c - 1) == 'W' and self._get_color(r, c + 1) == 'W')
                    is_vert_blocked = (self._get_color(r - 1, c) == 'B' and self._get_color(r + 1, c) == 'B')
                    if is_horiz_wall and is_vert_blocked:
                        return False
        return True

    def solve(self, k=0):
        """Recursively explores all 2^18 possible grid patterns."""
        if k == len(self.cells_to_fill):
            if self.check_word_lengths() and self.check_connectivity() and self.check_no_cheaters():
                self.valid_grids_count += 1
            return

        r, c = self.cells_to_fill[k]
        sr, sc = self.size - 1 - r, self.size - 1 - c

        # Option 1: Place black squares
        self.grid[r][c] = 'B'
        self.grid[sr][sc] = 'B'
        self.solve(k + 1)
        
        # Option 2: Place white squares
        self.grid[r][c] = 'W'
        self.grid[sr][sc] = 'W'
        self.solve(k + 1)
        
    def run(self):
        """Starts the solver and returns the final count."""
        self.solve()
        print(self.valid_grids_count)

if __name__ == '__main__':
    counter = CrosswordGridCounter(size=8)
    counter.run()
