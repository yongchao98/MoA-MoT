import sys

class CrosswordGridCounter:
    """
    This class counts the number of valid crossword grids for a given size,
    based on a specific set of rules:
    1. 180-degree rotational symmetry.
    2. All white squares are fully connected.
    3. Minimum word length (run of white squares) is 3.
    4. No 2x2 blocks of black squares.
    """

    def __init__(self, size):
        # The size must be even for this implementation of symmetry pairing.
        if size % 2 != 0:
            raise ValueError("Grid size must be even.")
        self.N = size
        self.grid = [[-1] * self.N for _ in range(self.N)]
        self.count = 0
        
        # We only need to decide the state of half the cells due to symmetry.
        # For an 8x8 grid, this is 32 cells. We can choose the top 4 rows.
        self.cells_to_fill = []
        for r in range(self.N // 2):
            for c in range(self.N):
                self.cells_to_fill.append((r, c))

    def solve(self):
        """
        Starts the backtracking search and prints the final count.
        Note: This is a computationally intensive search.
        """
        self.backtrack(0)
        # The problem asks for an equation, but this is a search problem.
        # The calculation is the search itself. The result is the final count.
        print(f"Number of squares to decide: {len(self.cells_to_fill)}")
        print(f"Total valid 8x8 crossword grids found: {self.count}")
        return self.count

    def backtrack(self, k):
        """
        Recursively explores grid configurations by filling one independent cell
        and its symmetric counterpart at each step.
        k: The index of the independent cell (from self.cells_to_fill) to decide.
        """
        if k == len(self.cells_to_fill):
            if self._is_valid_final_grid():
                self.count += 1
            return

        r, c = self.cells_to_fill[k]
        r_sym, c_sym = self.N - 1 - r, self.N - 1 - c

        # --- Option 1: Place a black square (1) ---
        self.grid[r][c] = 1
        self.grid[r_sym][c_sym] = 1
        
        # Pruning: If a 2x2 black block is formed, abandon this path.
        if not self._forms_2x2_black_square(r, c) and not self._forms_2x2_black_square(r_sym, c_sym):
            self.backtrack(k + 1)

        # --- Option 2: Place a white square (0) ---
        self.grid[r][c] = 0
        self.grid[r_sym][c_sym] = 0
        self.backtrack(k + 1)
        
        # Backtrack: reset the cell for the parent recursive call
        self.grid[r][c] = -1
        self.grid[r_sym][c_sym] = -1

    def _forms_2x2_black_square(self, r, c):
        """Checks if the black square at (r, c) completes any 2x2 block."""
        for dr_start in [-1, 0]:
            for dc_start in [-1, 0]:
                r_block, c_block = r + dr_start, c + dc_start
                if 0 <= r_block and r_block + 1 < self.N and 0 <= c_block and c_block + 1 < self.N:
                    if (self.grid[r_block][c_block] == 1 and
                        self.grid[r_block+1][c_block] == 1 and
                        self.grid[r_block][c_block+1] == 1 and
                        self.grid[r_block+1][c_block+1] == 1):
                        return True
        return False

    def _is_valid_final_grid(self):
        """Performs final validation checks on a fully generated grid."""
        white_squares = [(r, c) for r in range(self.N) for c in range(self.N) if self.grid[r][c] == 0]

        if not white_squares:
            return False

        # 1. Check for full interconnectivity of white squares
        q = [white_squares[0]]
        visited = {white_squares[0]}
        head = 0
        while head < len(q):
            r, c = q[head]
            head += 1
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < self.N and 0 <= nc < self.N and self.grid[nr][nc] == 0 and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))
        if len(visited) != len(white_squares):
            return False

        # 2. Check for minimum word length of 3
        for i in range(self.N):
            row_len = 0
            col_len = 0
            for j in range(self.N):
                # Check row
                if self.grid[i][j] == 0: row_len += 1
                else: 
                    if 0 < row_len < 3: return False
                    row_len = 0
                # Check column
                if self.grid[j][i] == 0: col_len += 1
                else:
                    if 0 < col_len < 3: return False
                    col_len = 0
            if 0 < row_len < 3 or 0 < col_len < 3:
                return False

        return True

if __name__ == '__main__':
    # Set a higher recursion limit for the deep search
    sys.setrecursionlimit(2000)
    counter = CrosswordGridCounter(8)
    counter.solve()