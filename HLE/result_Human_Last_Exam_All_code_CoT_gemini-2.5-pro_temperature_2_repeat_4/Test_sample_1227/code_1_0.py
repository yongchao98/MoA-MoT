import sys

# Set a higher recursion limit for the backtracking algorithm.
sys.setrecursionlimit(2000)

class CrosswordGridCounter:
    """
    Solves the crossword grid counting problem for a given size.
    """
    def __init__(self, size):
        if size % 2 != 0:
            raise ValueError("Grid size must be even for this symmetry implementation.")
        self.size = size
        self.grid = [[0] * size for _ in range(size)]
        self.solution_count = 0
        # We only need to iterate through the first half of the cells.
        self.unique_cells_count = (size * size) // 2

    def count_grids(self):
        """
        Starts the grid generation process and returns the final count.
        """
        self.generate(0)
        return self.solution_count

    def generate(self, index):
        """
        Recursively generates grid patterns using backtracking.
        """
        if index == self.unique_cells_count:
            if self.is_fully_valid():
                self.solution_count += 1
            return

        r, c = index // self.size, index % self.size
        sr, sc = self.size - 1 - r, self.size - 1 - c

        # Choice 1: Place a white square.
        # No pruning is needed for placing a white square.
        self.grid[r][c] = 0
        self.grid[sr][sc] = 0
        self.generate(index + 1)

        # Choice 2: Place a black square, with pruning.
        self.grid[r][c] = 1
        self.grid[sr][sc] = 1
        if self.is_placement_locally_valid(r, c):
            self.generate(index + 1)

    def is_placement_locally_valid(self, r, c):
        """
        Checks if placing a black square at (r, c) is valid based on
        already-placed squares (i.e., those with index smaller than current).
        This acts as a pruning mechanism.
        """
        # Pruning Rule 1: No 2x2 block of black squares ending at (r, c).
        if r > 0 and c > 0:
            if self.grid[r-1][c-1] == 1 and self.grid[r-1][c] == 1 and self.grid[r][c-1] == 1:
                return False

        # Pruning Rule 2: No short horizontal words created.
        # Checks for BWB or BWWB patterns to the left of the new black square.
        if c > 0 and self.grid[r][c - 1] == 0:  # If there's a white square to the left
            if c == 1:  # [Edge] W B
                return False
            if c > 1 and self.grid[r][c - 2] == 1:  # ...B W B
                return False
            if c == 2 and self.grid[r][c - 2] == 0:  # [Edge] W W B
                return False
            if c > 2 and self.grid[r][c - 3] == 1 and self.grid[r][c - 2] == 0:  # ...B W W B
                return False

        # Pruning Rule 3: No short vertical words created.
        if r > 0 and self.grid[r - 1][c] == 0:
            if r == 1:
                return False
            if r > 1 and self.grid[r - 2][c] == 1:
                return False
            if r == 2 and self.grid[r - 2][c] == 0:
                return False
            if r > 2 and self.grid[r - 3][c] == 1 and self.grid[r - 2][c] == 0:
                return False

        return True

    def is_fully_valid(self):
        """
        Performs a full validation of a completed grid.
        """
        # 1. Check for valid word lengths (>= 3).
        # Check rows
        for r in range(self.size):
            length = 0
            for c in range(self.size):
                if self.grid[r][c] == 0:
                    length += 1
                else:
                    if length in [1, 2]: return False
                    length = 0
            if length in [1, 2]: return False
        
        # Check columns
        for c in range(self.size):
            length = 0
            for r in range(self.size):
                if self.grid[r][c] == 0:
                    length += 1
                else:
                    if length in [1, 2]: return False
                    length = 0
            if length in [1, 2]: return False
            
        # 2. Check for connectivity of white squares.
        total_white = sum(row.count(0) for row in self.grid)
        if total_white == 0:
            return False

        q, visited = [], set()
        # Find first white square to start the search
        for r_start in range(self.size):
            try:
                c_start = self.grid[r_start].index(0)
                q.append((r_start, c_start))
                visited.add((r_start, c_start))
                break
            except ValueError:
                continue

        count_connected = 0
        while q:
            r, c = q.pop(0)
            count_connected += 1
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < self.size and 0 <= nc < self.size and \
                   self.grid[nr][nc] == 0 and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))

        if count_connected != total_white:
            return False

        # 3. Final check for "cheaters" (no 2x2 blocks of black squares).
        # This is partly redundant due to pruning, but necessary for correctness.
        for r in range(self.size - 1):
            for c in range(self.size - 1):
                if self.grid[r][c] == 1 and self.grid[r+1][c] == 1 and \
                   self.grid[r][c+1] == 1 and self.grid[r+1][c+1] == 1:
                    return False

        return True


if __name__ == '__main__':
    grid_size = 8
    counter = CrosswordGridCounter(grid_size)
    total_grids = counter.count_grids()
    print(total_grids)