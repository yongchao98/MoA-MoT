import collections

class CrosswordGridCounter:
    """
    This class calculates the number of valid 8x8 crossword grids based on a set of rules
    by using a pruned backtracking search algorithm.
    """

    def __init__(self, size=8):
        """Initializes the solver for a grid of a given size."""
        if size % 2 != 0:
            raise ValueError("Size must be even for this symmetry implementation.")
        self.size = size
        self.grid = [[-1] * size for _ in range(size)]
        self.count = 0
        # For an NxN grid with 180-degree symmetry, we only need to decide the
        # colors of the first N*N/2 cells. The other half is determined by symmetry.
        self.independent_cells = (size * size) // 2

    def solve(self):
        """
        Starts the recursive search for valid grids and returns the total count.
        """
        self.generate(0)
        return self.count

    def generate(self, k):
        """
        Recursively builds the grid by placing white (0) or black (1) squares.
        'k' is the index of the independent cell being decided (from 0 to 31).
        """
        # Base case: All 32 independent cells have been filled.
        if k == self.independent_cells:
            # Perform final checks that require a full grid.
            if self.check_final_constraints():
                self.count += 1
            return

        # Map k to grid coordinates (r, c) in the top half of the grid.
        r = k // self.size
        c = k % self.size
        # Determine the symmetric coordinates.
        sr, sc = self.size - 1 - r, self.size - 1 - c

        # Branch 1: Try placing a WHITE square.
        # Placing a white square cannot violate the min-word-length or no-cheater rules
        # at this local step, so we can always recurse.
        self.grid[r][c] = 0
        self.grid[sr][sc] = 0
        self.generate(k + 1)

        # Branch 2: Try placing a BLACK square.
        self.grid[r][c] = 1
        self.grid[sr][sc] = 1

        # After placing a black square, check for violations to prune this search branch.
        # We perform checks only when a full row is completed, as this is when violations
        # of the rules can be fully assessed.
        if c == self.size - 1:  # Check if row 'r' is now complete
            if self.is_row_pruning_valid(r):
                self.generate(k + 1)
        else: # If row is not yet complete, continue
            self.generate(k + 1)

    def is_row_pruning_valid(self, r):
        """Checks rules that can be evaluated once row 'r' and its symmetric part are complete."""
        sr = self.size - 1 - r

        # Check for words of length 1 or 2 in the completed row 'r' and its symmetric counterpart.
        if not self.check_line_word_length(self.grid[r]) or not self.check_line_word_length(self.grid[sr]):
            return False

        # Check for 2x2 "cheater" blocks. These can be checked when adjacent rows are complete.
        if r > 0:
            if not self.check_cheaters_in_row_pair(r - 1): # Checks between rows r-1 and r
                return False
        if sr < self.size - 1:
            if not self.check_cheaters_in_row_pair(sr): # Checks between rows sr and sr+1
                return False
        return True

    def check_final_constraints(self):
        """Performs checks that can only be done on a fully specified grid."""
        # Final cheater check for the middle rows, which is only possible now.
        if not self.check_cheaters_in_row_pair(self.size // 2 - 1):
            return False
        
        # Check all columns for valid word lengths.
        if not self.check_cols_word_length():
            return False
            
        # Check if all white squares are connected.
        if not self.check_connectivity():
            return False
            
        return True

    def check_line_word_length(self, line):
        """Checks a single list (row or column) for any sequence of 1 or 2 white squares."""
        length = 0
        for cell in line:
            if cell == 0:  # White square
                length += 1
            else:  # Black square
                if 0 < length < 3: # Found a word of length 1 or 2
                    return False
                length = 0
        if 0 < length < 3: # Check for a short word at the end of the line
            return False
        return True

    def check_cols_word_length(self):
        """Checks all grid columns for the minimum word length rule."""
        for c in range(self.size):
            col = [self.grid[r][c] for r in range(self.size)]
            if not self.check_line_word_length(col):
                return False
        return True

    def check_cheaters_in_row_pair(self, r):
        """Checks for any 2x2 blocks of black squares starting in row 'r'."""
        for c in range(self.size - 1):
            if (self.grid[r][c] == 1 and self.grid[r + 1][c] == 1 and
                    self.grid[r][c + 1] == 1 and self.grid[r + 1][c + 1] == 1):
                return False
        return True

    def check_connectivity(self):
        """Checks if all white squares form a single connected component using BFS."""
        white_squares = []
        for r in range(self.size):
            for c in range(self.size):
                if self.grid[r][c] == 0:
                    white_squares.append((r, c))

        # A valid crossword must have words, so an all-black grid is invalid.
        if not white_squares:
            return False
        
        # BFS to count reachable white squares from the first one found.
        q = collections.deque([white_squares[0]])
        visited = {white_squares[0]}
        
        count = 0
        while q:
            r_q, c_q = q.popleft()
            count += 1
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r_q + dr, c_q + dc
                if (0 <= nr < self.size and 0 <= nc < self.size and
                        self.grid[nr][nc] == 0 and (nr, nc) not in visited):
                    visited.add((nr, nc))
                    q.append((nr, nc))
        
        # If the number of visited squares equals the total, they are all connected.
        return count == len(white_squares)

if __name__ == '__main__':
    solver = CrosswordGridCounter(8)
    total_grids = solver.solve()
    print(total_grids)
    