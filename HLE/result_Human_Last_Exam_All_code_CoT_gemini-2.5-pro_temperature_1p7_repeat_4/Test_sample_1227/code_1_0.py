import sys

# We increase recursion limit for deep search paths, although pruning should prevent hitting it.
sys.setrecursionlimit(2000)

class CrosswordGridCounter:
    """
    Solves the crossword grid counting problem using a backtracking search.
    """

    def __init__(self, size):
        if size % 2 != 0:
            raise ValueError("Grid size must be even for this implementation.")
        self.N = size
        self.grid = [[0] * self.N for _ in range(self.N)]
        self.count = 0
        self.WHITE = 0
        self.BLACK = 1

    def solve(self):
        """Starts the backtracking process and returns the final count."""
        # We start filling from the top-left corner (0, 0).
        self.backtrack(0, 0)
        return self.count

    def check_line_word_lengths(self, line):
        """Checks a single list (a row or column) for word length violations."""
        length = 0
        for i in range(self.N):
            if line[i] == self.WHITE:
                length += 1
            else:
                if 1 <= length < 3:
                    return False
                length = 0
        if 1 <= length < 3: # Check for a word at the end of the line
            return False
        return True

    def check_all_column_lengths(self):
        """Checks all columns of the grid for word length violations."""
        for c in range(self.N):
            col = [self.grid[r][c] for r in range(self.N)]
            if not self.check_line_word_lengths(col):
                return False
        return True

    def check_connectivity(self):
        """Checks if all white squares are part of a single connected component."""
        white_squares = []
        for r in range(self.N):
            for c in range(self.N):
                if self.grid[r][c] == self.WHITE:
                    white_squares.append((r, c))

        # A valid puzzle must have words, so an all-black grid is invalid.
        if not white_squares:
            return False

        # Use BFS to find all reachable white squares from the first one.
        q = [white_squares[0]]
        visited = {white_squares[0]}
        head = 0
        while head < len(q):
            r, c = q[head]
            head += 1

            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < self.N and 0 <= nc < self.N and \
                   self.grid[nr][nc] == self.WHITE and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))
        
        # If the number of visited squares equals the total, they are all connected.
        return len(visited) == len(white_squares)

    def check_cheaters(self):
        """Checks for any 'cheater' black squares."""
        for r in range(self.N):
            for c in range(self.N):
                if self.grid[r][c] == self.BLACK:
                    has_top_white = r > 0 and self.grid[r-1][c] == self.WHITE
                    has_bottom_white = r < self.N - 1 and self.grid[r+1][c] == self.WHITE
                    is_walled_horizontally = (c == 0 or self.grid[r][c-1] == self.BLACK) and \
                                           (c == self.N - 1 or self.grid[r][c+1] == self.BLACK)
                    if has_top_white and has_bottom_white and is_walled_horizontally:
                        return False

                    has_left_white = c > 0 and self.grid[r][c-1] == self.WHITE
                    has_right_white = c < self.N - 1 and self.grid[r][c+1] == self.WHITE
                    is_walled_vertically = (r == 0 or self.grid[r-1][c] == self.BLACK) and \
                                           (r == self.N - 1 or self.grid[r+1][c] == self.BLACK)
                    if has_left_white and has_right_white and is_walled_vertically:
                        return False
        return True

    def backtrack(self, r, c):
        """
        Recursively fills the grid, pruning invalid branches.
        It only fills the top half (r < N/2), as the bottom is symmetric.
        """
        # Base Case: If we have filled the top half of the grid.
        if r == self.N // 2:
            # The grid is now fully determined by symmetry, so we run final validation.
            if self.check_all_column_lengths() and self.check_connectivity() and self.check_cheaters():
                self.count += 1
            return

        # Determine the next cell to fill.
        next_r, next_c = (r, c + 1) if c + 1 < self.N else (r + 1, 0)
        
        # --- Branch 1: Try placing a WHITE square ---
        self.grid[r][c] = self.WHITE
        self.grid[self.N - 1 - r][self.N - 1 - c] = self.WHITE
        
        # Pruning step: If a row is now complete, check its validity before recursing.
        if c == self.N - 1: 
            if self.check_line_word_lengths(self.grid[r]):
                self.backtrack(next_r, next_c)
        else: # Row is not yet complete, continue.
            self.backtrack(next_r, next_c)
            
        # --- Branch 2: Try placing a BLACK square ---
        self.grid[r][c] = self.BLACK
        self.grid[self.N - 1 - r][self.N - 1 - c] = self.BLACK

        # Pruning step: If a row is now complete, check its validity before recursing.
        if c == self.N - 1:
            if self.check_line_word_lengths(self.grid[r]):
                self.backtrack(next_r, next_c)
        else: # Row is not yet complete, continue.
            self.backtrack(next_r, next_c)

if __name__ == '__main__':
    grid_size = 8
    solver = CrosswordGridCounter(grid_size)
    num_grids = solver.solve()
    print(num_grids)