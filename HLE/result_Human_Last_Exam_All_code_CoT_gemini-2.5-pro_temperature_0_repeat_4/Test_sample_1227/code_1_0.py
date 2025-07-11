import collections

class CrosswordGridCounter:
    """
    Calculates the number of valid 8x8 crossword grids with a black border
    and 180-degree rotational symmetry.
    """
    def __init__(self, size=8):
        self.size = size
        self.inner_size = size - 2
        self.grid = [[-1] * size for _ in range(size)]
        self.count = 0
        
        # We only need to decide the colors for the first half of the inner grid cells.
        self.cells_to_fill = (self.inner_size * self.inner_size) // 2

        # Assume a standard black border. 1 represents a black square.
        for i in range(self.size):
            self.grid[0][i] = 1
            self.grid[self.size - 1][i] = 1
            self.grid[i][0] = 1
            self.grid[i][self.size - 1] = 1

    def solve(self):
        """
        Starts the recursive generation process and returns the final count.
        """
        self.generate(0)
        print(self.count)

    def generate(self, index):
        """
        Recursively fills the inner grid, exploring all possibilities.
        """
        # Base case: If all independent cells are filled, validate the grid.
        if index == self.cells_to_fill:
            if self.is_fully_valid():
                self.count += 1
            return

        # Map the linear index to 2D coordinates in the inner grid.
        r_inner = index // self.inner_size
        c_inner = index % self.inner_size
        
        # Map inner coordinates to the full grid coordinates.
        r = r_inner + 1
        c = c_inner + 1

        # Determine the symmetrically opposite coordinates.
        sr = self.size - 1 - r
        sc = self.size - 1 - c

        # --- Choice 1: Place a black square (1) ---
        self.grid[r][c] = 1
        self.grid[sr][sc] = 1
        self.generate(index + 1)

        # --- Choice 2: Place a white square (0) ---
        self.grid[r][c] = 0
        self.grid[sr][sc] = 0
        self.generate(index + 1)

        # Backtrack: Reset the cell for the parent recursive call.
        self.grid[r][c] = -1
        self.grid[sr][sc] = -1

    def is_fully_valid(self):
        """
        Checks if a fully generated grid meets all crossword rules.
        """
        return (self.check_no_2x2_black_blocks() and
                self.check_word_length() and
                self.check_connectivity())

    def check_no_2x2_black_blocks(self):
        """
        Ensures there are no 2x2 blocks of black squares.
        """
        for r in range(self.size - 1):
            for c in range(self.size - 1):
                if (self.grid[r][c] == 1 and self.grid[r+1][c] == 1 and
                    self.grid[r][c+1] == 1 and self.grid[r+1][c+1] == 1):
                    return False
        return True

    def check_word_length(self):
        """
        Ensures all horizontal and vertical words are at least 3 letters long.
        """
        # Check rows
        for r in range(self.size):
            length = 0
            for c in range(self.size):
                if self.grid[r][c] == 0:  # White square
                    length += 1
                else:  # Black square
                    if 0 < length < 3: return False
                    length = 0
            if 0 < length < 3: return False
        
        # Check columns
        for c in range(self.size):
            length = 0
            for r in range(self.size):
                if self.grid[r][c] == 0:  # White square
                    length += 1
                else:  # Black square
                    if 0 < length < 3: return False
                    length = 0
            if 0 < length < 3: return False
            
        return True

    def check_connectivity(self):
        """
        Ensures all white squares are connected using Breadth-First Search (BFS).
        """
        white_squares = []
        for r in range(self.size):
            for c in range(self.size):
                if self.grid[r][c] == 0:
                    white_squares.append((r, c))

        if not white_squares:
            return True  # An all-black grid is trivially valid.

        q = collections.deque([white_squares[0]])
        visited = {white_squares[0]}
        
        while q:
            r, c = q.popleft()
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if (0 <= nr < self.size and 0 <= nc < self.size and
                    self.grid[nr][nc] == 0 and (nr, nc) not in visited):
                    visited.add((nr, nc))
                    q.append((nr, nc))

        return len(visited) == len(white_squares)

if __name__ == '__main__':
    counter = CrosswordGridCounter(size=8)
    counter.solve()