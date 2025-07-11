import collections
import sys

class CrosswordCounter:
    """
    A class to count the number of valid crossword grids of a given size.
    It uses a backtracking algorithm to explore all symmetric patterns and
    validates them against a set of crossword rules.
    """
    def __init__(self, size=8):
        # A grid of size N is represented as (N+2)x(N+2) to have a permanent
        # black border (1=black, 0=white), which simplifies boundary checks.
        self.size = size
        self.grid = [[1] * (self.size + 2) for _ in range(self.size + 2)]
        self.count = 0
        self.independent_cells = []
        self._prepare_grid_and_cells()

    def _prepare_grid_and_cells(self):
        """
        Initializes the grid with a black border and identifies the set of
        independent cells to decide during the search, based on symmetry.
        """
        # A standard crossword convention is a solid black outer border.
        # We enforce this for the 8x8 puzzle grid itself.
        for i in range(self.size):
            self.grid[1][i + 1] = 1
            self.grid[self.size][i + 1] = 1
            self.grid[i + 1][1] = 1
            self.grid[i + 1][self.size] = 1

        # Due to 180-degree symmetry, we only need to decide the color for
        # one cell in each symmetric pair. We find these "independent" cells.
        visited = set()
        # Iterate over the inner grid (from (2,2) to (7,7) in the padded 10x10 grid)
        for r in range(2, self.size):
            for c in range(2, self.size):
                if (r, c) not in visited:
                    self.independent_cells.append((r, c))
                    sym_r = (self.size + 1) - r
                    sym_c = (self.size + 1) - c
                    visited.add((r, c))
                    visited.add((sym_r, sym_c))

    def _is_connected(self):
        """Checks if all white squares in the grid are connected."""
        white_squares = [(r, c) for r in range(1, self.size + 1) for c in range(1, self.size + 1) if self.grid[r][c] == 0]

        if not white_squares:
            return False  # No white squares means no words.

        q = collections.deque([white_squares[0]])
        visited = {white_squares[0]}
        
        while q:
            r, c = q.popleft()
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if self.grid[nr][nc] == 0 and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))
        
        return len(visited) == len(white_squares)

    def _get_word_info(self, temp_grid):
        """
        Calculates the number of words in a grid and checks for min length validity.
        A word is a maximal contiguous block of white squares.
        """
        word_count = 0
        is_valid = True
        # Horizontal scan
        for r in range(1, self.size + 1):
            c = 1
            while c <= self.size:
                if temp_grid[r][c] == 0 and temp_grid[r][c - 1] == 1:
                    word_len = sum(1 for k in range(c, self.size + 1) if temp_grid[r][k] == 0)
                    word_count += 1
                    if word_len < 3: is_valid = False
                    c += word_len
                else: c += 1
        
        # Vertical scan
        for c in range(1, self.size + 1):
            r = 1
            while r <= self.size:
                if temp_grid[r][c] == 0 and temp_grid[r - 1][c] == 1:
                    word_len = sum(1 for k in range(r, self.size + 1) if temp_grid[k][c] == 0)
                    word_count += 1
                    if is_valid and word_len < 3: is_valid = False
                    r += word_len
                else: r += 1
        return word_count, is_valid

    def _has_no_cheaters(self):
        """
        Checks if the grid contains any 'cheater' squares. A cheater is a
        black square that, when flipped to white, results in a valid grid
        with the same number of words.
        """
        black_squares = [(r, c) for r in range(2, self.size) for c in range(2, self.size) if self.grid[r][c] == 1]
        
        if not black_squares: return True

        original_word_count, _ = self._get_word_info(self.grid)

        for r_b, c_b in black_squares:
            temp_grid = [row[:] for row in self.grid]
            temp_grid[r_b][c_b] = 0
            
            flipped_word_count, is_flipped_valid = self._get_word_info(temp_grid)
            
            if is_flipped_valid and flipped_word_count == original_word_count:
                return False  # Found a cheater square

        return True

    def _solve_recursive(self, k=0):
        """
        Recursively explores patterns by setting independent cells and their
        symmetric pairs to black or white.
        """
        if k == len(self.independent_cells):
            # Base case: full pattern is generated. Now, validate it.
            _, is_len_valid = self._get_word_info(self.grid)
            if is_len_valid and self._is_connected() and self._has_no_cheaters():
                self.count += 1
            return

        r, c = self.independent_cells[k]
        sym_r = (self.size + 1) - r
        sym_c = (self.size + 1) - c

        # Branch 1: Try placing white squares
        self.grid[r][c], self.grid[sym_r][sym_c] = 0, 0
        self._solve_recursive(k + 1)

        # Branch 2: Try placing black squares
        self.grid[r][c], self.grid[sym_r][sym_c] = 1, 1
        self._solve_recursive(k + 1)

    def run(self):
        """Runs the solver and prints the result."""
        print("Calculating the number of valid 8x8 crossword grids...", file=sys.stderr)
        print("This involves checking 2^18 = 262,144 base patterns.", file=sys.stderr)
        print("This may take a minute...", file=sys.stderr)
        
        self._solve_recursive()
        
        final_count = self.count
        print(f"The final number of possible grids is: {final_count}")
        return final_count

if __name__ == '__main__':
    counter = CrosswordCounter(size=8)
    counter.run()