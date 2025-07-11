import sys

# It is recommended to run this script using PyPy for a significant speed improvement.
# If using CPython, it might take a few minutes.

class CrosswordGridCounter:
    """
    A class to count valid crossword grids based on a set of rules.
    """
    def __init__(self, size):
        if size % 2 != 0:
            raise ValueError("Grid size must be even for this symmetry implementation.")
        self.size = size
        self.grid = [[-1] * size for _ in range(size)]
        self.valid_grid_count = 0
        self.cells_to_decide = self._get_cells_to_decide()

    def _get_cells_to_decide(self):
        """
        Generates the list of cells in the inner grid that need a color decision.
        Due to symmetry, we only need to decide for half of the cells.
        For an 8x8 grid with a black border, this is the top half of the inner 6x6 grid.
        Rows 1, 2, 3 and Columns 1 through 6. (18 cells)
        """
        cells = []
        inner_start = 1
        inner_end = self.size - 1
        # Iterate through the top half of the inner grid's rows
        for r in range(inner_start, self.size // 2 + 1):
            # Iterate through all columns of the inner grid
            for c in range(inner_start, inner_end):
                # For the middle row pair, only iterate through half the columns
                if r == self.size // 2 and c >= self.size // 2:
                    continue
                cells.append((r, c))
        return cells

    def solve(self):
        """
        Starts the process of finding and counting all valid grids.
        """
        # Set outer border to black, a standard crossword convention.
        for i in range(self.size):
            self.grid[0][i] = 1
            self.grid[self.size - 1][i] = 1
            self.grid[i][0] = 1
            self.grid[i][self.size - 1] = 1
        
        self.backtrack(0)
        return self.valid_grid_count

    def backtrack(self, k):
        """
        Recursively explores all possible grid patterns.
        """
        if k == len(self.cells_to_decide):
            if self._is_grid_valid():
                self.valid_grid_count += 1
            return

        r, c = self.cells_to_decide[k]
        r_sym = self.size - 1 - r
        c_sym = self.size - 1 - c

        # Try placing black squares
        self.grid[r][c] = 1
        self.grid[r_sym][c_sym] = 1
        self.backtrack(k + 1)

        # Try placing white squares
        self.grid[r][c] = 0
        self.grid[r_sym][c_sym] = 0
        self.backtrack(k + 1)

    def _is_grid_valid(self):
        """
        Checks a completed grid against all rules.
        """
        # The order of checks is chosen to fail fast on cheaper tests.
        if not self._has_no_short_words():
            return False
        if self._contains_cheater_squares():
            return False
        if not self._is_fully_connected():
            return False
        return True

    def _has_no_short_words(self):
        """
        Checks for any words (runs of white squares) of length 1 or 2.
        """
        for i in range(self.size):
            h_len = 0
            v_len = 0
            for j in range(self.size):
                # Check horizontal word length
                if self.grid[i][j] == 0:
                    h_len += 1
                else:
                    if 0 < h_len < 3: return False
                    h_len = 0
                # Check vertical word length
                if self.grid[j][i] == 0:
                    v_len += 1
                else:
                    if 0 < v_len < 3: return False
                    v_len = 0
            # Check at the end of a row/col
            if 0 < h_len < 3: return False
            if 0 < v_len < 3: return False
        return True

    def _is_fully_connected(self):
        """
        Checks if all white squares form a single contiguous region using BFS.
        """
        white_squares = []
        for r in range(self.size):
            for c in range(self.size):
                if self.grid[r][c] == 0:
                    white_squares.append((r, c))

        if not white_squares:
            return True

        q = [white_squares[0]]
        visited = {white_squares[0]}
        head = 0
        while head < len(q):
            r, c = q[head]
            head += 1
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < self.size and 0 <= nc < self.size and \
                   self.grid[nr][nc] == 0 and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))
        
        return len(visited) == len(white_squares)

    def _contains_cheater_squares(self):
        """
        A grid contains a cheater if a black square's removal would not change the word count.
        This occurs if flipping it to white doesn't merge any existing words.
        """
        for r in range(1, self.size - 1):
            for c in range(1, self.size - 1):
                if self.grid[r][c] == 1:  # If it's a black square
                    # Check if it would merge horizontal words if flipped
                    merges_h = (self.grid[r][c-1] == 0 and self.grid[r][c+1] == 0)
                    # Check if it would merge vertical words if flipped
                    merges_v = (self.grid[r-1][c] == 0 and self.grid[r+1][c] == 0)
                    
                    if not merges_h and not merges_v:
                        return True # Found a cheater square
        return False


if __name__ == '__main__':
    grid_size = 8
    counter = CrosswordGridCounter(grid_size)
    total_grids = counter.solve()
    print(f"For a {grid_size}x{grid_size} puzzle with the specified rules (180-degree symmetry, min word length 3, full interconnect, no cheaters, and black borders), there are {total_grids} possible grids.")
    print("<<<" + str(total_grids) + ">>>")