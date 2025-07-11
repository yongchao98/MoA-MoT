import sys

class CrosswordGridCounter:
    """
    Solves the problem of counting valid 8x8 crossword grids based on a set of rules.
    """
    def __init__(self, size):
        self.N = size
        # Grid state: -1 for undecided, 0 for white, 1 for black
        self.grid = [[-1] * self.N for _ in range(self.N)]
        self.count = 0
        # Determine the cells to fill, leveraging symmetry.
        # We only need to iterate through the first half of the cells.
        self.cells_to_fill = []
        for r in range(self.N):
            for c in range(self.N):
                if (r, c) <= (self.N - 1 - r, self.N - 1 - c):
                    self.cells_to_fill.append((r, c))

    def run(self):
        """Starts the solving process and prints the result."""
        self._solve_recursive(0)
        # The puzzle asks for the final answer. This code calculates it.
        # Final result is printed to the console.
        print(f"Number of possible 8x8 crossword grids: {self.count}")
    
    def _get_color(self, r, c):
        """Gets the color of a square, treating the border as black."""
        if not (0 <= r < self.N and 0 <= c < self.N):
            return 1  # Border is considered black
        return self.grid[r][c]

    def _solve_recursive(self, k):
        """The main backtracking function."""
        # Base case: if all independent cells are filled, validate the full grid.
        if k == len(self.cells_to_fill):
            if self._is_grid_fully_valid():
                self.count += 1
            return

        r, c = self.cells_to_fill[k]
        sr, sc = self.N - 1 - r, self.N - 1 - c

        # Branch 1: Try placing WHITE squares
        self.grid[r][c] = 0
        if (r, c) != (sr, sc): self.grid[sr][sc] = 0
        self._solve_recursive(k + 1)

        # Branch 2: Try placing BLACK squares
        self.grid[r][c] = 1
        if (r, c) != (sr, sc): self.grid[sr][sc] = 1
        
        # Prune this branch if the placement is invalid
        if self._is_placement_locally_valid(r, c) and ((r, c) == (sr, sc) or self._is_placement_locally_valid(sr, sc)):
            self._solve_recursive(k + 1)

        # Backtrack: reset the cells for other branches of the search
        self.grid[r][c] = -1
        if (r, c) != (sr, sc): self.grid[sr][sc] = -1

    def _is_placement_locally_valid(self, r, c):
        """
        Pruning function: checks if placing a black square at (r, c) creates
        an immediate, unrecoverable violation of the rules.
        """
        # Rule: No cheater squares.
        # Check if the new square at (r,c) is itself a cheater.
        if self._get_color(r, c - 1) == 1 and self._get_color(r, c + 1) == 1: return False
        if self._get_color(r - 1, c) == 1 and self._get_color(r + 1, c) == 1: return False
        
        # Check if the new square makes an existing neighbor black square a cheater.
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc
            if self._get_color(nr, nc) == 1:
                # Check horizontal neighbors of the neighbor
                if self._get_color(nr, nc - 1) == 1 and self._get_color(nr, nc + 1) == 1: return False
                # Check vertical neighbors of the neighbor
                if self._get_color(nr - 1, nc) == 1 and self._get_color(nr + 1, nc) == 1: return False
        
        # Rule: Minimum word length of 3.
        # Check if we just enclosed a short run of white squares.
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            run_len = 0
            rr, cc = r + dr, c + dc
            while self._get_color(rr, cc) == 0:
                run_len += 1
                rr += dr
                cc += dc
            # If the run is now bordered by black squares and is too short, prune.
            if self._get_color(rr, cc) == 1 and 0 < run_len < 3:
                return False
                
        return True

    def _is_grid_fully_valid(self):
        """
        Performs all final validation checks on a completed grid.
        """
        # 1. Connectivity check
        white_squares = [(r, c) for r in range(self.N) for c in range(self.N) if self.grid[r][c] == 0]
        if not white_squares:
            return False  # A valid puzzle must have some words.
        
        q = [white_squares[0]]
        visited = {white_squares[0]}
        head = 0
        while head < len(q):
            r, c = q[head]; head += 1
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < self.N and 0 <= nc < self.N and self.grid[nr][nc] == 0 and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))
        if len(visited) != len(white_squares):
            return False

        # 2. Minimum word length and cheater checks
        for r in range(self.N):
            for c in range(self.N):
                # Check for short horizontal words starting at (r, c)
                if self._get_color(r, c - 1) == 1 and self._get_color(r, c) == 0:
                    run = 0
                    while self._get_color(r, c + run) == 0: run += 1
                    if self._get_color(r, c + run) == 1 and run < 3: return False
                # Check for short vertical words starting at (r, c)
                if self._get_color(r - 1, c) == 1 and self._get_color(r, c) == 0:
                    run = 0
                    while self._get_color(r + run, c) == 0: run += 1
                    if self._get_color(r + run, c) == 1 and run < 3: return False
                # Check for cheater squares
                if self.grid[r][c] == 1:
                    if self._get_color(r, c - 1) == 1 and self._get_color(r, c + 1) == 1: return False
                    if self._get_color(r - 1, c) == 1 and self._get_color(r + 1, c) == 1: return False
        
        return True

if __name__ == '__main__':
    # Increasing Python's recursion limit may be necessary for larger grid sizes,
    # but should be fine for 8x8.
    sys.setrecursionlimit(2000)
    counter = CrosswordGridCounter(size=8)
    counter.run()
