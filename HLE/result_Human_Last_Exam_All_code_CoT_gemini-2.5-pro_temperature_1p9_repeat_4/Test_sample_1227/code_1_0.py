class CrosswordCounter:
    """
    Calculates the number of valid 8x8 crossword grids based on a set of rules.

    The rules are:
    1. 8x8 grid size.
    2. 180-degree rotational symmetry.
    3. Minimum word length of 3 (no runs of 1 or 2 white squares).
    4. All white squares must be fully interconnected.
    5. No "cheater" squares, interpreted as no 2x2 blocks of black squares.
    """
    def __init__(self, size=8):
        self.size = size
        self.grid = [[0] * size for _ in range(size)] # 0 for white, 1 for black
        self.count = 0
        
        # We only need to decide the color for the top half of the grid (32 squares for 8x8).
        # The bottom half is determined by the 180-degree symmetry rule.
        self.coords_to_fill = []
        for r in range(self.size // 2):
            for c in range(self.size):
                self.coords_to_fill.append((r, c))

    def solve(self):
        """Public method to start the solving process."""
        self._fill(0)
        return self.count

    def _fill(self, k):
        """Recursively fills the grid, trying both black and white for each square."""
        # Base Case: All unique squares have been filled. The grid is ready for final validation.
        if k == len(self.coords_to_fill):
            if self._is_final_grid_valid():
                self.count += 1
            return

        r, c = self.coords_to_fill[k]
        rs, cs = self.size - 1 - r, self.size - 1 - c

        # --- Choice 1: Place a BLACK square ---
        self.grid[r][c] = 1
        self.grid[rs][cs] = 1
        
        # Pruning: If this move creates a 2x2 black block, abandon this path.
        if not self._check_2x2_violation_at(r, c) and not self._check_2x2_violation_at(rs, cs):
            self._fill(k + 1)

        # Backtrack: Revert the change to explore the other choice.
        self.grid[r][c] = 0
        self.grid[rs][cs] = 0

        # --- Choice 2: Place a WHITE square ---
        # The grid square is already white, so we just proceed to the next step.
        self._fill(k + 1)

    def _check_2x2_violation_at(self, r, c):
        """Checks if the newly placed black square at (r, c) completes a 2x2 black block."""
        # A square at (r,c) can be part of four potential 2x2 blocks (as a corner).
        for dr in [-1, 0]:
            for dc in [-1, 0]:
                r_start, c_start = r + dr, c + dc
                if 0 <= r_start < self.size - 1 and 0 <= c_start < self.size - 1:
                    if all(self.grid[nr][nc] == 1 
                           for nr in [r_start, r_start + 1] 
                           for nc in [c_start, c_start + 1]):
                        return True
        return False

    def _is_final_grid_valid(self):
        """Performs expensive validation checks on a fully generated grid."""
        # Check for trivial grids (all black or all white).
        total_black = sum(sum(row) for row in self.grid)
        if total_black == 0 or total_black == self.size * self.size:
            return False

        if not self._is_connected(total_black):
            return False
            
        if not self._has_valid_word_lengths():
            return False
            
        return True

    def _is_connected(self, total_black):
        """Checks if all white squares are connected using BFS."""
        total_white = self.size * self.size - total_black
        
        q = []
        visited = set()
        
        # Find the first white square to start the search.
        for r_start in range(self.size):
            for c_start in range(self.size):
                if self.grid[r_start][c_start] == 0:
                    q.append((r_start, c_start))
                    visited.add((r_start, c_start))
                    break
            if q:
                break
        
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
        
        return len(visited) == total_white
        
    def _has_valid_word_lengths(self):
        """Checks if all horizontal and vertical words are at least 3 letters long."""
        # Check rows
        for r in range(self.size):
            length = 0
            for c in range(self.size):
                if self.grid[r][c] == 0:
                    length += 1
                else:
                    if 0 < length < 3: return False
                    length = 0
            if 0 < length < 3: return False
        
        # Check columns
        for c in range(self.size):
            length = 0
            for r in range(self.size):
                if self.grid[r][c] == 0:
                    length += 1
                else:
                    if 0 < length < 3: return False
                    length = 0
            if 0 < length < 3: return False
            
        return True

if __name__ == '__main__':
    # This calculation is computationally intensive and may take a few minutes to run.
    solver = CrosswordCounter(size=8)
    result = solver.solve()
    print(result)