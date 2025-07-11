import collections

class CrosswordGridCounter:
    """
    This class calculates the number of valid 8x8 crossword grids with
    180-degree rotational symmetry based on a set of common puzzle constraints.
    """
    def __init__(self, size=8):
        if size % 2 != 0:
            raise ValueError("Size must be even for this symmetry model.")
        self.size = size
        self.grid = [[-1] * size for _ in range(size)]
        self.count = 0

        # Assumption: The grid is enclosed in a border of black squares.
        # This simplifies the problem to filling the inner (size-2)x(size-2) grid.
        for i in range(self.size):
            self.grid[0][i] = self.grid[self.size - 1][i] = 1
            self.grid[i][0] = self.grid[i][self.size - 1] = 1

        # Define the independent squares we need to decide.
        # For an 8x8 grid with a border, this is the top half of the inner 6x6 grid.
        self.independent_squares = []
        inner_size = self.size - 2
        for r_offset in range(inner_size // 2):
            for c_offset in range(inner_size):
                self.independent_squares.append((r_offset + 1, c_offset + 1))

    def solve(self):
        """Starts the recursive search and returns the total count of valid grids."""
        self._recursive_solve(0)
        return self.count

    def _recursive_solve(self, k):
        """
        Recursively fills the grid by choosing black or white for each independent square.
        """
        # Base case: if all independent squares are filled, we have a complete grid.
        if k == len(self.independent_squares):
            # Now, validate the complete grid against the rules.
            if self._check_cheaters() and self._check_word_length() and self._check_connectivity():
                self.count += 1
            return

        r, c = self.independent_squares[k]
        r_sym = self.size - 1 - r
        c_sym = self.size - 1 - c

        # Branch 1: Try placing a black square.
        self.grid[r][c] = 1
        self.grid[r_sym][c_sym] = 1
        self._recursive_solve(k + 1)

        # Branch 2: Try placing a white square.
        self.grid[r][c] = 0
        self.grid[r_sym][c_sym] = 0
        self._recursive_solve(k + 1)

    def _check_cheaters(self):
        """Checks for 'cheater' squares, proxied by 'no 2x2 black blocks'."""
        for r in range(self.size - 1):
            for c in range(self.size - 1):
                if (self.grid[r][c] == 1 and
                    self.grid[r + 1][c] == 1 and
                    self.grid[r][c + 1] == 1 and
                    self.grid[r + 1][c + 1] == 1):
                    return False
        return True

    def _check_line_word_length(self, line):
        """Helper to check word lengths in a single row or column."""
        length = 0
        for square in line:
            if square == 0:  # White square
                length += 1
            else:  # Black square
                if 0 < length < 3:
                    return False
                length = 0
        if 0 < length < 3:  # Check for a word at the end of the line
            return False
        return True

    def _check_word_length(self):
        """Checks that all words (horizontal/vertical) are at least 3 letters long."""
        # Check all rows
        for r in range(self.size):
            if not self._check_line_word_length(self.grid[r]):
                return False
        # Check all columns
        for c in range(self.size):
            col = [self.grid[r][c] for r in range(self.size)]
            if not self._check_line_word_length(col):
                return False
        return True

    def _check_connectivity(self):
        """Checks if all white squares are connected into a single component."""
        white_squares = []
        for r in range(self.size):
            for c in range(self.size):
                if self.grid[r][c] == 0:
                    white_squares.append((r, c))

        if not white_squares:
            return False  # A grid with no white squares is not a valid puzzle.

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
    total_grids = counter.solve()
    print(f"The number of possible 8x8 crossword grids is: {total_grids}")