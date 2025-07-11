import sys

class CrosswordGridCounter:
    """
    This class calculates the number of possible 8x8 crossword grids
    based on a specific set of rules:
    1.  8x8 grid size.
    2.  180-degree rotational symmetry.
    3.  All words (contiguous white squares) must be at least 3 letters long.
    4.  All white squares must be fully interconnected.
    5.  The grid is bordered by black squares (an assumed standard convention).
    6.  No 2x2 blocks of black squares (a common, concrete interpretation of 'no cheaters').
    """

    def __init__(self, size=8):
        if size % 2 != 0:
            raise ValueError("Size must be even for this symmetry implementation.")
        self.size = size
        # -1: undecided, 0: white, 1: black
        self.grid = [[-1] * size for _ in range(size)]
        self.count = 0
        
        # Assumption: The grid border is made of black squares.
        for i in range(size):
            self.grid[0][i] = 1
            self.grid[size - 1][i] = 1
            self.grid[i][0] = 1
            self.grid[i][size - 1] = 1

        # List of squares in the inner grid that we need to decide.
        # Due to symmetry, we only need to decide for the first half.
        self.squares_to_decide = []
        inner_size = size - 2
        for r_inner in range(inner_size):
            for c_inner in range(inner_size):
                # An inner square (r_inner, c_inner) corresponds to (r,c) in the full grid
                r, c = r_inner + 1, c_inner + 1
                
                # To pick one from each symmetric pair, we only take the one
                # that comes first in a row-by-row scan of the grid.
                if r * self.size + c <= (self.size - 1 - r) * self.size + (self.size - 1 - c):
                    self.squares_to_decide.append((r, c))

    def solve(self):
        """Starts the recursive grid generation process and returns the total count."""
        self._generate(0)
        return self.count

    def _generate(self, k):
        """
        Recursively tries placing white and black squares, backtracking through
        the possibilities.
        """
        # Base case: All undecided squares have been filled.
        if k == len(self.squares_to_decide):
            if self._is_valid_grid():
                self.count += 1
            return

        r, c = self.squares_to_decide[k]
        r_sym, c_sym = self.size - 1 - r, self.size - 1 - c

        # --- Option 1: Place a WHITE square ---
        self.grid[r][c] = 0
        if r != r_sym or c != c_sym: # Avoid double-writing for center squares
            self.grid[r_sym][c_sym] = 0
        self._generate(k + 1)

        # --- Option 2: Place a BLACK square ---
        self.grid[r][c] = 1
        if r != r_sym or c != c_sym:
            self.grid[r_sym][c_sym] = 1
            
        # Pruning: Check if this placement creates a 2x2 black block.
        if not self._has_2x2_black_block_at(r, c) and not self._has_2x2_black_block_at(r_sym, c_sym):
            self._generate(k + 1)
        
        # Backtracking is implicit as we overwrite the grid in subsequent calls.

    def _is_valid_grid(self):
        """Checks if a fully generated grid meets all the criteria."""
        if not self._check_word_lengths():
            return False
        if not self._check_connectivity():
            return False
        # The 2x2 block check is mostly handled by pruning, but a final check is safest.
        if self._has_any_2x2_black_block():
            return False
        return True

    def _check_word_lengths(self):
        # Check rows for words shorter than 3
        for r in range(self.size):
            length = 0
            for c in range(self.size):
                if self.grid[r][c] == 0:  # White
                    length += 1
                else:  # Black
                    if 0 < length < 3: return False
                    length = 0
            if 0 < length < 3: return False

        # Check columns for words shorter than 3
        for c in range(self.size):
            length = 0
            for r in range(self.size):
                if self.grid[r][c] == 0:
                    length += 1
                else:
                    if 0 < length < 3: return False
                    length = 0
            if 0 < length < 3: return False
        return True

    def _check_connectivity(self):
        total_white = sum(row.count(0) for row in self.grid)
        if total_white == 0:
            return True 

        q = []
        visited = set()
        
        first_white = next(((r, c) for r in range(self.size) for c in range(self.size) if self.grid[r][c] == 0), None)
        
        if not first_white: return True

        q.append(first_white)
        visited.add(first_white)
        count_connected = 0
        
        while q:
            r, c = q.pop(0)
            count_connected += 1
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < self.size and 0 <= nc < self.size and self.grid[nr][nc] == 0 and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))

        return total_white == count_connected

    def _has_2x2_black_block_at(self, r, c):
        """Checks if placing a black square at (r, c) completes a 2x2 black block."""
        for dr in [-1, 0]:
            for dc in [-1, 0]:
                r_start, c_start = r + dr, c + dc
                if r_start >= 0 and c_start >= 0 and r_start + 1 < self.size and c_start + 1 < self.size:
                    if (self.grid[r_start][c_start] == 1 and
                        self.grid[r_start + 1][c_start] == 1 and
                        self.grid[r_start][c_start + 1] == 1 and
                        self.grid[r_start + 1][c_start + 1] == 1):
                        return True
        return False
        
    def _has_any_2x2_black_block(self):
        """Checks the entire grid for any 2x2 black blocks."""
        for r in range(self.size - 1):
            for c in range(self.size - 1):
                if (self.grid[r][c] == 1 and self.grid[r+1][c] == 1 and
                    self.grid[r][c+1] == 1 and self.grid[r+1][c+1] == 1):
                    return True
        return False

counter = CrosswordGridCounter(size=8)
total_grids = counter.solve()
print(total_grids)