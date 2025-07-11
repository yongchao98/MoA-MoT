import sys

class CrosswordGridCounter:
    """
    A class to count valid crossword grids based on a set of rules.
    """
    def __init__(self, size):
        if size % 2 != 0:
            raise ValueError("Grid size must be even for this implementation.")
        self.size = size
        self.grid = [[-1] * size for _ in range(size)]
        self.count = 0
        # We only need to decide the color for the first half of the cells
        # due to 180-degree rotational symmetry.
        self.cells_to_fill = []
        for r in range(size // 2):
            for c in range(size):
                self.cells_to_fill.append((r, c))

    def count_grids(self):
        """
        Starts the recursive process to count all valid grids.
        """
        self._generate(0)
        return self.count

    def _is_valid(self):
        """
        Checks if a fully generated grid satisfies all crossword rules.
        """
        # Rule: No "cheater" squares (interpreted as no 2x2 black squares)
        for r in range(self.size - 1):
            for c in range(self.size - 1):
                if (self.grid[r][c] == 1 and
                    self.grid[r+1][c] == 1 and
                    self.grid[r][c+1] == 1 and
                    self.grid[r+1][c+1] == 1):
                    return False

        # Rule: Minimum word length of 3
        # Check horizontal words
        for r in range(self.size):
            for c in range(self.size):
                if self.grid[r][c] == 0: # If it's a white square
                    # Check if it's the start of a horizontal word
                    if c == 0 or self.grid[r][c-1] == 1:
                        length = 0
                        k = c
                        while k < self.size and self.grid[r][k] == 0:
                            length += 1
                            k += 1
                        if length < 3:
                            return False
        # Check vertical words
        for c in range(self.size):
            for r in range(self.size):
                if self.grid[r][c] == 0: # If it's a white square
                    # Check if it's the start of a vertical word
                    if r == 0 or self.grid[r-1][c] == 1:
                        length = 0
                        k = r
                        while k < self.size and self.grid[k][c] == 0:
                            length += 1
                            k += 1
                        if length < 3:
                            return False

        # Rule: Full interconnectivity of white squares
        total_white = sum(row.count(0) for row in self.grid)
        if total_white == 0:
            return True # A grid with no white squares is valid by these rules

        # Find the first white square to start the search from
        start_node = None
        for r in range(self.size):
            for c in range(self.size):
                if self.grid[r][c] == 0:
                    start_node = (r, c)
                    break
            if start_node:
                break
        
        # Use BFS to find all connected white squares
        q = [start_node]
        visited = {start_node}
        connected_count = 0
        while q:
            r, c = q.pop(0)
            connected_count += 1
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if (0 <= nr < self.size and 0 <= nc < self.size and
                        self.grid[nr][nc] == 0 and (nr, nc) not in visited):
                    visited.add((nr, nc))
                    q.append((nr, nc))

        return connected_count == total_white

    def _generate(self, k):
        """
        Recursively fills the grid and checks for valid solutions.
        """
        # Base case: if all unique cells are filled, validate the full grid
        if k == len(self.cells_to_fill):
            if self._is_valid():
                self.count += 1
            return

        r, c = self.cells_to_fill[k]
        r_sym, c_sym = self.size - 1 - r, self.size - 1 - c

        # Try placing a white square (0)
        self.grid[r][c] = 0
        self.grid[r_sym][c_sym] = 0
        self._generate(k + 1)

        # Try placing a black square (1)
        self.grid[r][c] = 1
        self.grid[r_sym][c_sym] = 1
        self._generate(k + 1)

if __name__ == '__main__':
    # Note: This computation is very intensive and will take a long time to complete.
    grid_size = 8
    counter = CrosswordGridCounter(grid_size)
    total_grids = counter.count_grids()
    print(total_grids)