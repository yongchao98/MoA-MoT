import sys

class CrosswordGridCounter:
    """
    This class calculates the number of valid crossword grids of a given size
    based on a set of rules including 180-degree rotational symmetry.
    """

    def __init__(self, size):
        """
        Initializes the counter. We assume the grid is bounded by black squares,
        so we only need to fill the inner (size-2)x(size-2) grid.
        """
        if size % 2 != 0 or size < 4:
            raise ValueError("Size must be an even integer of at least 4.")
        
        self.size = size
        # Start with a grid that has a 1-square thick black border.
        self.grid = [[1] * self.size for _ in range(self.size)]
        self.count = 0
        
        inner_size = self.size - 2
        
        # Determine the list of independent cells in the inner grid to iterate over.
        # For an 8x8 grid, the inner grid is 6x6, and we fill its top half (3 rows).
        self.cells_to_fill = []
        for r_inner in range(inner_size // 2):
            for c_inner in range(inner_size):
                # Convert inner coordinates to full grid coordinates
                self.cells_to_fill.append((r_inner + 1, c_inner + 1))
        
        # For a 6x6 inner grid, there are (6*6)/2 = 18 independent cells.
        self.num_independent_cells = len(self.cells_to_fill)

    def solve(self):
        """
        Starts the recursive generation and validation process.
        """
        self.generate(0)
        return self.count

    def generate(self, k):
        """
        Recursively generates all possible symmetric patterns for the inner grid.
        k is the index of the independent cell we are currently deciding.
        """
        # Base case: If we have filled all independent cells, the grid is complete.
        if k == self.num_independent_cells:
            # A complete symmetric grid is formed, now we validate it.
            if self.is_valid_grid():
                self.count += 1
            return

        # Get the coordinates of the current independent cell
        r, c = self.cells_to_fill[k]
        # Get the coordinates of its symmetric counterpart
        r_sym, c_sym = self.size - 1 - r, self.size - 1 - c

        # Recursive step: Try both black and white for the cell pair.

        # Option 1: Place white squares (0)
        self.grid[r][c] = 0
        self.grid[r_sym][c_sym] = 0
        self.generate(k + 1)

        # Option 2: Place black squares (1)
        self.grid[r][c] = 1
        self.grid[r_sym][c_sym] = 1
        self.generate(k + 1)

    def is_valid_grid(self):
        """
        Checks if a fully generated grid meets all crossword puzzle criteria.
        """
        # --- Check 1: Word Lengths (must be >= 3) ---
        # Horizontal words
        for r in range(self.size):
            c = 0
            while c < self.size:
                if self.grid[r][c] == 0:
                    length = 0
                    start_c = c
                    while c < self.size and self.grid[r][c] == 0:
                        length += 1
                        c += 1
                    if length < 3:
                        return False
                else:
                    c += 1
        
        # Vertical words
        for c in range(self.size):
            r = 0
            while r < self.size:
                if self.grid[r][c] == 0:
                    length = 0
                    start_r = r
                    while r < self.size and self.grid[r][c] == 0:
                        length += 1
                        r += 1
                    if length < 3:
                        return False
                else:
                    r += 1

        # --- Check 2: Connectivity of White Squares ---
        total_white_squares = sum(row.count(0) for row in self.grid)
        if total_white_squares == 0:
            return True # A grid of all black squares is valid.

        q = []
        visited = set()
        # Find the first white square to start the search from.
        for r_start in range(self.size):
            try:
                c_start = self.grid[r_start].index(0)
                q.append((r_start, c_start))
                visited.add((r_start, c_start))
                break
            except ValueError:
                continue
        
        # Perform a BFS to find all connected white squares.
        head = 0
        while head < len(q):
            curr_r, curr_c = q[head]
            head += 1
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                next_r, next_c = curr_r + dr, curr_c + dc
                if 0 <= next_r < self.size and 0 <= next_c < self.size and \
                   self.grid[next_r][next_c] == 0 and (next_r, next_c) not in visited:
                    visited.add((next_r, next_c))
                    q.append((next_r, next_c))
        
        if len(visited) != total_white_squares:
            return False

        # --- Check 3: No "Cheater" Squares ---
        # A black square is a "cheater" if it's surrounded on all 4 sides
        # by other black squares or the grid boundary.
        for r in range(self.size):
            for c in range(self.size):
                if self.grid[r][c] == 1:
                    up_blocked = (r == 0) or (self.grid[r-1][c] == 1)
                    down_blocked = (r == self.size - 1) or (self.grid[r+1][c] == 1)
                    left_blocked = (c == 0) or (self.grid[r][c-1] == 1)
                    right_blocked = (c == self.size - 1) or (self.grid[r][c+1] == 1)
                    
                    if up_blocked and down_blocked and left_blocked and right_blocked:
                        return False

        return True

if __name__ == '__main__':
    # For an 8x8 puzzle, we run the counter with size=8.
    counter = CrosswordGridCounter(8)
    result = counter.solve()
    print(f"The number of possible 8x8 crossword grids is: {result}")