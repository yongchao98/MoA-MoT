import sys

# Increase recursion limit for deep search paths.
sys.setrecursionlimit(2000)

class CrosswordGridCounter:
    """
    This class systematically generates and validates 8x8 crossword grids
    based on a set of rules, including rotational symmetry, minimum word length,
    connectivity, and the absence of "cheater" squares.
    """
    def __init__(self, size=8):
        self.N = size
        self.grid = [[-1] * self.N for _ in range(self.N)]
        self.solution_count = 0
        self.half_size = (self.N * self.N) // 2

    def solve(self):
        """
        Starts the recursive generation process and prints the final count.
        """
        self.generate(0)
        print("Final calculation:")
        print(f"Number of valid 8x8 crossword grids = {self.solution_count}")
    
    def is_in_bounds(self, r, c):
        """Checks if coordinates are within the grid."""
        return 0 <= r < self.N and 0 <= c < self.N

    def generate(self, index):
        """
        Recursively places black or white squares, pruning invalid branches.
        """
        if index == self.half_size:
            # A full symmetric grid pattern has been generated. Now, validate it.
            if self.is_fully_valid():
                self.solution_count += 1
            return

        r, c = index // self.N, index % self.N

        # --- Try placing a black square (1) ---
        self.grid[r][c] = 1
        self.grid[self.N - 1 - r][self.N - 1 - c] = 1
        
        # Pruning: Check if this placement creates an immediate invalid word.
        # This check is crucial for performance.
        can_be_black = True
        # Check horizontal segment to the left, which is now terminated.
        if c > 0 and self.grid[r][c - 1] == 0:
            length = sum(1 for i in range(c - 1, -1, -1) if self.grid[r][i] == 0)
            if 0 < length < 3:
                can_be_black = False
        
        # Check vertical segment above, which is now terminated.
        if can_be_black and r > 0 and self.grid[r - 1][c] == 0:
            length = sum(1 for i in range(r - 1, -1, -1) if self.grid[i][c] == 0)
            if 0 < length < 3:
                can_be_black = False

        if can_be_black:
            self.generate(index + 1)

        # --- Try placing a white square (0) ---
        self.grid[r][c] = 0
        self.grid[self.N - 1 - r][self.N - 1 - c] = 0
        self.generate(index + 1)

    def is_fully_valid(self):
        """
        Runs all validation checks on a completed grid.
        """
        if not self.check_final_word_lengths(): return False
        if not self.check_connectivity(): return False
        if not self.check_for_cheaters(): return False
        return True

    def check_final_word_lengths(self):
        """Ensures no words of length 1 or 2 exist anywhere in the grid."""
        # Check all rows
        for r in range(self.N):
            length = 0
            for c in range(self.N):
                if self.grid[r][c] == 0: length += 1
                else:
                    if 0 < length < 3: return False
                    length = 0
            if 0 < length < 3: return False
        
        # Check all columns
        for c in range(self.N):
            length = 0
            for r in range(self.N):
                if self.grid[r][c] == 0: length += 1
                else:
                    if 0 < length < 3: return False
                    length = 0
            if 0 < length < 3: return False
        
        return True

    def check_connectivity(self):
        """Ensures all white squares are connected."""
        white_squares = {(r, c) for r in range(self.N) for c in range(self.N) if self.grid[r][c] == 0}
        
        # An all-black grid is invalid. An all-white grid is also invalid (no words).
        if not white_squares or len(white_squares) == self.N * self.N:
            return False

        q = [next(iter(white_squares))]
        visited = {q[0]}
        
        head = 0
        while head < len(q):
            r, c = q[head]
            head += 1
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if self.is_in_bounds(nr, nc) and self.grid[nr][nc] == 0 and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))
                    
        return len(visited) == len(white_squares)

    def _is_row_ok(self, r):
        """Helper to check a single row's word lengths."""
        length = 0
        for c in range(self.N):
            if self.grid[r][c] == 0: length += 1
            else:
                if 0 < length < 3: return False
                length = 0
        return not (0 < length < 3)

    def _is_col_ok(self, c):
        """Helper to check a single column's word lengths."""
        length = 0
        for r in range(self.N):
            if self.grid[r][c] == 0: length += 1
            else:
                if 0 < length < 3: return False
                length = 0
        return not (0 < length < 3)

    def check_for_cheaters(self):
        """
        A black square is a cheater if flipping it to white (a) results
        in a valid grid (re word length) and (b) does not change the word count.
        """
        for r in range(self.N):
            for c in range(self.N):
                if self.grid[r][c] == 1: # For each black square
                    # Check conditions for merging words
                    merges_h = self.is_in_bounds(r, c-1) and self.grid[r][c-1] == 0 and \
                               self.is_in_bounds(r, c+1) and self.grid[r][c+1] == 0
                    merges_v = self.is_in_bounds(r-1, c) and self.grid[r-1][c] == 0 and \
                               self.is_in_bounds(r+1, c) and self.grid[r+1][c] == 0
                    
                    # Tentatively flip the square
                    self.grid[r][c] = 0
                    
                    # A flip is valid only if it doesn't create invalid words in its row/col
                    is_valid_flip = self._is_row_ok(r) and self._is_col_ok(c)
                    
                    # Flip back to restore grid state
                    self.grid[r][c] = 1

                    if is_valid_flip and not merges_h and not merges_v:
                        return False # Found a cheater square, grid is invalid
        return True

if __name__ == '__main__':
    counter = CrosswordGridCounter(size=8)
    counter.solve()
<<<116>>>