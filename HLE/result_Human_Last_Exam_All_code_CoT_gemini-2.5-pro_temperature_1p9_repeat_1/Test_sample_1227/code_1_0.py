import sys

class CrosswordGridCounter:
    """
    Calculates the number of valid 8x8 crossword grids based on a set of rules.
    
    The core logic assumes that the grid is bounded by black squares, a standard
    convention in American-style crosswords. This reduces the search space to a
    manageable size.
    """
    def __init__(self, size=8):
        if size % 2 != 0:
            raise ValueError("Grid size must be even for this symmetric logic.")
        self.size = size
        self.grid = [[1] * size for _ in range(size)] # Initialize with all black squares
        self.count = 0
        # The coordinates of the inner cells in the top-left quadrant that we need to decide.
        self.cells_to_set = []
        inner_size = size - 2
        for i in range(inner_size // 2):
            for j in range(inner_size):
                self.cells_to_set.append((i + 1, j + 1))
        
    def _is_connected(self):
        """Checks if all white squares form a single connected component using BFS."""
        white_squares = []
        for r in range(self.size):
            for c in range(self.size):
                if self.grid[r][c] == 0:
                    white_squares.append((r, c))

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
                if 0 <= nr < self.size and 0 <= nc < self.size and \
                   self.grid[nr][nc] == 0 and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))
        return len(visited) == len(white_squares)

    def _has_valid_word_lengths(self):
        """Checks if all word lengths are 3 or more."""
        # Check rows
        for r in range(self.size):
            length = 0
            for c in range(self.size):
                if self.grid[r][c] == 0:
                    length += 1
                else:
                    if 1 <= length <= 2: return False
                    length = 0
            if 1 <= length <= 2: return False
        
        # Check columns
        for c in range(self.size):
            length = 0
            for r in range(self.size):
                if self.grid[r][c] == 0:
                    length += 1
                else:
                    if 1 <= length <= 2: return False
                    length = 0
            if 1 <= length <= 2: return False
            
        return True

    def _has_no_cheaters(self):
        """Checks that every black square has at least one white neighbor."""
        for r in range(self.size):
            for c in range(self.size):
                if self.grid[r][c] == 1: # Black square
                    has_white_neighbor = False
                    for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                        nr, nc = r + dr, c + dc
                        if 0 <= nr < self.size and 0 <= nc < self.size and self.grid[nr][nc] == 0:
                            has_white_neighbor = True
                            break
                    if not has_white_neighbor:
                        return False
        return True

    def _is_valid_grid(self):
        """Runs all validation checks."""
        return self._has_valid_word_lengths() and \
               self._has_no_cheaters() and \
               self._is_connected()

    def solve(self, k=0):
        """Recursively explores all possible grid configurations."""
        if k == len(self.cells_to_set):
            if self._is_valid_grid():
                self.count += 1
            return

        r, c = self.cells_to_set[k]
        
        # Trial 1: Place a black square
        # The grid is already black, so we don't need to set it, just proceed.
        self.grid[r][c] = 1
        self.grid[self.size - 1 - r][self.size - 1 - c] = 1
        self.solve(k + 1)
        
        # Trial 2: Place a white square
        self.grid[r][c] = 0
        self.grid[self.size - 1 - r][self.size - 1 - c] = 0
        self.solve(k + 1)
        
        # Backtrack: Restore grid to its all-black state for the parent call
        self.grid[r][c] = 1
        self.grid[self.size - 1 - r][self.size - 1 - c] = 1

    def run(self):
        """Starts the calculation and prints the result."""
        self.solve()
        print(f"Total number of possible 8x8 crossword grids: {self.count}")
        # The problem asks to output the number from the final equation.
        # As there is no equation, we output the only relevant number.
        print("Final Answer:")
        print(self.count)


if __name__ == '__main__':
    counter = CrosswordGridCounter(size=8)
    counter.run()