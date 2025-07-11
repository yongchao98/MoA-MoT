import sys

class CrosswordGridCounter:
    """
    A class to count valid 8x8 crossword grids based on a set of rules.

    The rules are:
    1. 180-degree rotational symmetry.
    2. The outermost border is black.
    3. All words (sequences of white squares) must be at least 3 letters long.
    4. All white squares must be fully interconnected.
    5. No 2x2 blocks of black squares (as a proxy for the "no cheaters" rule).
    """

    def __init__(self, size=8):
        """
        Initializes the counter with a grid of a given size.
        """
        if size % 2 != 0:
            raise ValueError("Size must be an even number for this symmetry implementation.")
        
        self.size = size
        self.grid = [[None for _ in range(size)] for _ in range(size)]
        self.valid_grid_count = 0

        # We only need to decide the colors for the inner grid.
        self.inner_size = size - 2
        
        # Due to symmetry, we only need to decide on half of the inner grid cells.
        self.num_cells_to_decide = (self.inner_size * self.inner_size) // 2
        
        # Pre-calculate the coordinates of the cells we need to decide.
        # We choose the top half of the inner grid.
        self.cells_to_decide_coords = []
        for i in range(self.num_cells_to_decide):
            inner_r = i // self.inner_size
            inner_c = i % self.inner_size
            # Map inner grid coordinates to the full grid coordinates (with a 1-square border)
            self.cells_to_decide_coords.append((inner_r + 1, inner_c + 1))

    def count_grids(self):
        """
        Starts the grid generation and counting process.
        Returns the total count of valid grids.
        """
        # Rule: The outer border must be black (1 = black, 0 = white).
        for i in range(self.size):
            self.grid[i][0] = 1
            self.grid[i][self.size - 1] = 1
            self.grid[0][i] = 1
            self.grid[self.size - 1][i] = 1
            
        self._generate_recursively(0)
        return self.valid_grid_count

    def _generate_recursively(self, k):
        """
        Recursively tries all patterns for the undecided cells using backtracking.
        k: The index of the current cell we are placing.
        """
        # Base case: All independent cells have been set.
        if k == self.num_cells_to_decide:
            # The grid is fully populated, so we validate it.
            if self._is_valid_grid():
                self.valid_grid_count += 1
            return

        # Get coordinates for the current cell and its symmetric counterpart.
        r, c = self.cells_to_decide_coords[k]
        sym_r, sym_c = self.size - 1 - r, self.size - 1 - c

        # --- Recursive Step ---
        
        # Branch 1: Try placing white squares (0).
        self.grid[r][c] = 0
        self.grid[sym_r][sym_c] = 0
        self._generate_recursively(k + 1)
        
        # Branch 2: Try placing black squares (1).
        self.grid[r][c] = 1
        self.grid[sym_r][sym_c] = 1
        self._generate_recursively(k + 1)

    def _is_valid_grid(self):
        """
        Checks a fully generated grid against all the puzzle rules.
        """
        # First, find all white squares. If none, it's not a valid puzzle.
        white_squares = []
        for r in range(self.size):
            for c in range(self.size):
                if self.grid[r][c] == 0:
                    white_squares.append((r, c))
        
        if not white_squares:
            return False

        # Rule: Full interconnectivity of white squares.
        if not self._check_connectivity(white_squares):
            return False

        # Rule: Minimum word length of 3.
        if not self._check_word_lengths():
            return False
            
        # Rule: No 2x2 blocks of black squares.
        if not self._check_no_2x2_black_squares():
            return False
            
        return True

    def _check_connectivity(self, white_squares):
        """
        Checks if all white squares are connected using Breadth-First Search (BFS).
        """
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

    def _check_word_lengths(self):
        """
        Checks that all horizontal and vertical words have a length of at least 3.
        """
        # Check rows for horizontal words
        for r in range(self.size):
            length = 0
            for c in range(self.size):
                if self.grid[r][c] == 0: # Part of a potential word
                    length += 1
                else:  # Black square marks the end of a word
                    if 0 < length < 3:
                        return False
                    length = 0
            if 0 < length < 3: # Check for a word ending at the grid edge
                return False

        # Check columns for vertical words
        for c in range(self.size):
            length = 0
            for r in range(self.size):
                if self.grid[r][c] == 0:
                    length += 1
                else:  # Black square
                    if 0 < length < 3:
                        return False
                    length = 0
            if 0 < length < 3: # Check for a word ending at the grid edge
                return False
        return True
    
    def _check_no_2x2_black_squares(self):
        """
        Checks that there are no 2x2 blocks of black squares.
        """
        for r in range(self.size - 1):
            for c in range(self.size - 1):
                if self.grid[r][c] == 1 and \
                   self.grid[r+1][c] == 1 and \
                   self.grid[r][c+1] == 1 and \
                   self.grid[r+1][c+1] == 1:
                    return False
        return True

if __name__ == '__main__':
    # Increase recursion limit to be safe, though default is likely sufficient.
    sys.setrecursionlimit(2000)
    
    counter = CrosswordGridCounter(size=8)
    total_valid_grids = counter.count_grids()
    print(total_valid_grids)