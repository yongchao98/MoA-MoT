import sys

class CrosswordGridCounter:
    """
    This class calculates the number of valid 8x8 crossword grids based on a set of rules.
    """
    def __init__(self, size):
        self.N = size
        if self.N % 2 != 0:
            # The symmetry logic is simpler for even-sized grids.
            raise ValueError("This implementation is designed for even-sized grids.")
        self.grid = [[-1] * self.N for _ in range(self.N)] # -1 represents an unfilled square
        self.count = 0

    def solve(self):
        """
        Starts the backtracking search and returns the total count of valid grids.
        """
        # We start filling the grid from the first pair of squares.
        # k=0 represents the square at (0,0) and its symmetric partner.
        self.backtrack(0)
        return self.count

    def check_word_length(self, arr):
        """
        Checks a single row or column for the minimum word length of 3.
        A 'word' is a contiguous sequence of white squares (0).
        A black square is represented by 1.
        """
        line_str = ''.join(map(str, arr))
        # Split the line by black squares ('1') to find the words.
        words = line_str.split('1')
        for word in words:
            # Any sequence of 1 or 2 white squares is invalid.
            if 1 <= len(word) <= 2:
                return False
        return True

    def check_no_2x2_black_squares(self, grid):
        """Checks for any 2x2 blocks of black squares (1) in the final grid."""
        for r in range(self.N - 1):
            for c in range(self.N - 1):
                if grid[r][c] == 1 and grid[r+1][c] == 1 and \
                   grid[r][c+1] == 1 and grid[r+1][c+1] == 1:
                    return False
        return True

    def check_connectivity(self, grid):
        """
        Checks if all white squares (0) are connected in a single component using BFS.
        """
        white_squares = []
        for r in range(self.N):
            for c in range(self.N):
                if grid[r][c] == 0:
                    white_squares.append((r, c))

        if not white_squares:
            # An all-black grid has no white squares to be disconnected.
            # It passes this test, but will be filtered by the no_2x2_black_squares check.
            return True
        
        # Start a Breadth-First Search (BFS) from the first white square found.
        q = [white_squares[0]]
        visited = {white_squares[0]}
        head = 0
        
        while head < len(q):
            r, c = q[head]
            head += 1
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < self.N and 0 <= nc < self.N and \
                   grid[nr][nc] == 0 and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))
        
        # If the number of visited squares equals the total number of white squares, they are all connected.
        return len(visited) == len(white_squares)

    def backtrack(self, k):
        """
        Recursively fills the grid, prunes invalid branches, and counts valid solutions.
        k: the index of the square-pair to fill, from 0 to (N*N/2 - 1).
        """
        # Base case: The grid is fully populated.
        if k == (self.N * self.N) // 2:
            # Row word lengths were checked during pruning. Now check columns.
            for c in range(self.N):
                col = [self.grid[r][c] for r in range(self.N)]
                if not self.check_word_length(col):
                    return

            # Check for connectivity and 2x2 black squares.
            if self.check_connectivity(self.grid) and self.check_no_2x2_black_squares(self.grid):
                self.count += 1
            return

        # Map k to grid coordinates (r, c) in the upper half of the grid.
        r, c = k // self.N, k % self.N
        # Get coordinates of the symmetrically opposite square.
        r_sym, c_sym = self.N - 1 - r, self.N - 1 - c

        # ---- Branch 1: Try placing black squares (1) ----
        self.grid[r][c] = 1
        self.grid[r_sym][c_sym] = 1
        
        prune_branch = False
        # Pruning check: If a row has just been completed, check its word length.
        if (k + 1) % self.N == 0:
            if not self.check_word_length(self.grid[r]):
                 prune_branch = True

        if not prune_branch:
            self.backtrack(k + 1)

        # ---- Branch 2: Try placing white squares (0) ----
        self.grid[r][c] = 0
        self.grid[r_sym][c_sym] = 0

        prune_branch = False
        if (k + 1) % self.N == 0:
            if not self.check_word_length(self.grid[r]):
                 prune_branch = True

        if not prune_branch:
            self.backtrack(k + 1)
        
        # Restore state for parent call in recursion tree.
        self.grid[r][c] = -1
        self.grid[r_sym][c_sym] = -1


if __name__ == '__main__':
    grid_size = 8
    counter = CrosswordGridCounter(grid_size)
    num_grids = counter.solve()
    print(num_grids)