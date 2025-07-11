import sys

class CrosswordGridCounter:
    """
    This class finds the number of valid 8x8 crossword grids with 180-degree
    rotational symmetry based on a set of standard rules.
    """
    def __init__(self, size=8):
        if size % 2 != 0 or size < 4:
            raise ValueError("Size must be an even integer of at least 4.")
        
        self.size = size
        # We will use 1 for black squares and 0 for white squares.
        self.grid = [[-1] * size for _ in range(size)]
        self.count = 0
        self.cells_to_fill = []

        # Assumption 1: The grid has black borders.
        for i in range(size):
            self.grid[i][0] = 1
            self.grid[i][size - 1] = 1
            self.grid[0][i] = 1
            self.grid[size - 1][i] = 1
            
        # Due to symmetry, we only need to decide the colors for the cells
        # in the upper half of the inner grid. For an 8x8 grid, these are
        # the cells (r, c) where 1 <= r <= 3 and 1 <= c <= 6.
        for r in range(1, size // 2):
            for c in range(1, size - 1):
                self.cells_to_fill.append((r, c))

    def solve(self):
        """
        Starts the generation and counting process.
        """
        self._generate(0)
        # The problem asks how many possible grids can be made.
        print(f"Number of possible 8x8 crossword grids: {self.count}")

    def _generate(self, k):
        """
        Recursively generates grid patterns by filling one unique cell at a time.
        """
        if k == len(self.cells_to_fill):
            if self._is_valid_grid():
                self.count += 1
            return

        r, c = self.cells_to_fill[k]
        
        # Option 1: Place a white square
        self.grid[r][c] = 0
        self.grid[self.size - 1 - r][self.size - 1 - c] = 0
        self._generate(k + 1)

        # Option 2: Place a black square
        self.grid[r][c] = 1
        self.grid[self.size - 1 - r][self.size - 1 - c] = 1
        self._generate(k + 1)

    def _is_valid_grid(self):
        """
        Checks if a fully generated grid is valid according to all rules.
        """
        # Rule 1: No "cheater" squares (interpreted as no 2x2 black blocks)
        for r in range(self.size - 1):
            for c in range(self.size - 1):
                if (self.grid[r][c] == 1 and self.grid[r+1][c] == 1 and
                    self.grid[r][c+1] == 1 and self.grid[r+1][c+1] == 1):
                    return False

        # Rule 2: Minimum word length of 3
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

        # Rule 3: All white squares must be fully interconnected
        white_squares = []
        for r in range(self.size):
            for c in range(self.size):
                if self.grid[r][c] == 0:
                    white_squares.append((r,c))

        if not white_squares:
            return False # An all-black grid isn't a puzzle.

        # BFS to check for connectivity
        q = [white_squares[0]]
        visited = {white_squares[0]}
        head = 0
        while head < len(q):
            r_curr, c_curr = q[head]
            head += 1
            for dr, dc in [(0,1), (0,-1), (1,0), (-1,0)]:
                nr, nc = r_curr + dr, c_curr + dc
                # No need to check bounds since borders are black (and thus not in white_squares)
                if self.grid[nr][nc] == 0 and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))

        if len(visited) != len(white_squares):
            return False
            
        return True

if __name__ == '__main__':
    # Increasing the recursion limit is a good practice for deep recursion,
    # though it might not be strictly necessary for this problem size.
    sys.setrecursionlimit(2000)
    
    counter = CrosswordGridCounter(size=8)
    counter.solve()