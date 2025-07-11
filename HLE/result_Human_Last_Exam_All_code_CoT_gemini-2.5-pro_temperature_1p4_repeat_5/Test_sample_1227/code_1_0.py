import sys

class CrosswordGridCounter:
    """
    This class calculates the number of valid crossword grids of a given size
    that satisfy a specific set of rules:
    1. 180-degree rotational symmetry.
    2. Minimum word length of 3 (for all horizontal and vertical words).
    3. Full interconnectivity of all white squares.
    4. No 2x2 blocks of black squares (a common interpretation of "no cheaters").
    """

    def __init__(self, size):
        if size % 2 != 0:
            raise ValueError("Grid size must be even for this symmetry implementation.")
        self.size = size
        # Grid state: 0 for white, 1 for black, -1 for undecided
        self.grid = [[-1] * size for _ in range(size)]
        self.count = 0

    def solve(self):
        """
        Starts the backtracking process to find all valid grid patterns.
        The search space is the set of 2^(N*N/2) symmetric grids.
        """
        # index from 0 to (size*size/2 - 1) corresponds to the unique cells
        # in the top half of the grid that determine the entire pattern.
        self.backtrack(0)
        return self.count

    def backtrack(self, index):
        """
        Recursively explores placing white or black squares, pruning invalid branches.
        """
        # Base case: if we have filled the top half of the grid.
        if index == (self.size * self.size) // 2:
            # The grid is fully determined. Perform final checks.
            # Row checks have been done via pruning. Now check columns.
            if not self._check_all_column_word_lengths():
                return
            # Finally, check for connectivity of white squares.
            if self._check_connectivity():
                self.count += 1
            return

        r = index // self.size
        c = index % self.size
        r_sym = self.size - 1 - r
        c_sym = self.size - 1 - c

        # --- Choice 1: Place White Squares ---
        self.grid[r][c] = 0
        self.grid[r_sym][c_sym] = 0
        # If this placement completes a row, we can check word length rule for that row.
        can_recurse_white = True
        if c == self.size - 1:
            if not (self._check_row_word_length(r) and self._check_row_word_length(r_sym)):
                can_recurse_white = False
        if can_recurse_white:
            self.backtrack(index + 1)

        # --- Choice 2: Place Black Squares ---
        self.grid[r][c] = 1
        self.grid[r_sym][c_sym] = 1
        # Pruning: Check if this placement creates a 2x2 black square violation.
        can_recurse_black = True
        if self._violates_2x2_rule_at(r, c) or self._violates_2x2_rule_at(r_sym, c_sym):
            can_recurse_black = False
        else:
            # If valid so far and it completes a row, check word length rule.
            if c == self.size - 1:
                if not (self._check_row_word_length(r) and self._check_row_word_length(r_sym)):
                    can_recurse_black = False
        
        if can_recurse_black:
            self.backtrack(index + 1)

        # Backtrack: reset the squares to undecided for other search branches.
        self.grid[r][c] = -1
        self.grid[r_sym][c_sym] = -1

    def _violates_2x2_rule_at(self, r, c):
        """Checks if placing a black square at (r, c) creates a 2x2 block."""
        # A new 2x2 block must involve (r,c). Check all 4 possible 2x2 blocks.
        # Note: self.grid can have -1 (undecided), so we must check for == 1.
        # Check block where (r, c) is the bottom-right corner
        if r > 0 and c > 0 and self.grid[r-1][c] == 1 and self.grid[r][c-1] == 1 and self.grid[r-1][c-1] == 1:
            return True
        # Check block where (r, c) is the bottom-left corner
        if r > 0 and c < self.size - 1 and self.grid[r-1][c] == 1 and self.grid[r][c+1] == 1 and self.grid[r-1][c+1] == 1:
            return True
        # Check block where (r, c) is the top-right corner
        if r < self.size - 1 and c > 0 and self.grid[r+1][c] == 1 and self.grid[r][c-1] == 1 and self.grid[r+1][c-1] == 1:
            return True
        # Check block where (r, c) is the top-left corner
        if r < self.size - 1 and c < self.size - 1 and self.grid[r+1][c] == 1 and self.grid[r][c+1] == 1 and self.grid[r+1][c+1] == 1:
            return True
        return False

    def _check_row_word_length(self, r):
        """Checks if all horizontal words in a given row are of length >= 3."""
        length = 0
        for c in range(self.size):
            if self.grid[r][c] == 0:  # White square
                length += 1
            else:  # Black square
                if 0 < length < 3: return False
                length = 0
        if 0 < length < 3: return False  # Check at the end of the row
        return True

    def _check_all_column_word_lengths(self):
        """Checks if all vertical words in the grid are of length >= 3."""
        for c in range(self.size):
            length = 0
            for r in range(self.size):
                if self.grid[r][c] == 0:  # White square
                    length += 1
                else:  # Black square
                    if 0 < length < 3: return False
                    length = 0
            if 0 < length < 3: return False  # Check at the end of the column
        return True

    def _check_connectivity(self):
        """Checks if all white squares form a single connected component using BFS."""
        white_squares = []
        for r in range(self.size):
            for c in range(self.size):
                if self.grid[r][c] == 0:
                    white_squares.append((r, c))
        
        # A valid crossword must have at least one white square.
        if not white_squares:
            return False

        q = [white_squares[0]]
        visited = {white_squares[0]}
        head = 0
        while head < len(q):
            r, c = q[head]
            head += 1
            
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if (0 <= nr < self.size and 0 <= nc < self.size and
                        self.grid[nr][nc] == 0 and (nr, nc) not in visited):
                    visited.add((nr, nc))
                    q.append((nr, nc))
                    
        return len(visited) == len(white_squares)


# The main execution part of the script
if __name__ == '__main__':
    size = 8
    counter = CrosswordGridCounter(size)
    result = counter.solve()
    print(result)