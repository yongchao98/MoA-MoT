import collections

class CrosswordGridCounter:
    """
    A class to count valid crossword grids based on a set of rules.
    """
    def __init__(self, size):
        if size % 2 != 0:
            raise ValueError("Grid size must be even for this implementation.")
        self.N = size
        self.grid = [[-1] * self.N for _ in range(self.N)]
        self.count = 0

    def count_grids(self):
        """
        Starts the recursive grid generation process and returns the total count.
        """
        # We only need to decide the first N*N/2 squares due to symmetry.
        self._generate(0)
        return self.count

    def _generate(self, index):
        """
        Recursively generates grid patterns by placing white/black squares.
        """
        # Base case: If the first half of the grid is filled, validate the full grid.
        if index >= (self.N * self.N) // 2:
            if self._is_fully_valid():
                self.count += 1
            return

        r = index // self.N
        c = index % self.N
        sym_r, sym_c = self.N - 1 - r, self.N - 1 - c

        # --- Option 1: Place white squares ---
        self.grid[r][c] = 0
        self.grid[sym_r][sym_c] = 0
        self._generate(index + 1)

        # --- Option 2: Place black squares ---
        self.grid[r][c] = 1
        self.grid[sym_r][sym_c] = 1
        # Pruning step: check for 2x2 black blocks immediately.
        if self._is_safe_to_place_black(r, c) and self._is_safe_to_place_black(sym_r, sym_c):
            self._generate(index + 1)

    def _is_safe_to_place_black(self, r, c):
        """
        Checks if placing a black square at (r, c) creates a 2x2 block of black squares.
        Assumes grid[r][c] is already 1.
        """
        # Check the 4 possible 2x2 squares that (r,c) can be a corner of.
        for dr in [-1, 0]:
            for dc in [-1, 0]:
                is_block = True
                for i in range(2):
                    for j in range(2):
                        nr, nc = r + dr + i, c + dc + j
                        if not (0 <= nr < self.N and 0 <= nc < self.N and self.grid[nr][nc] == 1):
                            is_block = False
                            break
                    if not is_block:
                        break
                if is_block:
                    return False
        return True

    def _is_fully_valid(self):
        """
        Checks if a fully generated grid meets all constraints.
        """
        return self._check_connectivity() and self._check_word_lengths()

    def _check_connectivity(self):
        """
        Checks if all white squares are connected using a Breadth-First Search (BFS).
        """
        total_white = sum(row.count(0) for row in self.grid)
        if total_white == 0:
            return False

        q = collections.deque()
        visited = set()
        
        # Find the first white square to start the search
        for r_start in range(self.N):
            for c_start in range(self.N):
                if self.grid[r_start][c_start] == 0:
                    q.append((r_start, c_start))
                    visited.add((r_start, c_start))
                    break
            if q:
                break
        
        count_reachable = 0
        while q:
            r, c = q.popleft()
            count_reachable += 1
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < self.N and 0 <= nc < self.N and \
                   self.grid[nr][nc] == 0 and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))

        return count_reachable == total_white

    def _check_word_lengths(self):
        """
        Checks that all horizontal and vertical words have a minimum length of 3.
        """
        for i in range(self.N):
            h_len, v_len = 0, 0
            for j in range(self.N):
                # Horizontal check
                if self.grid[i][j] == 0:
                    h_len += 1
                else:
                    if 0 < h_len < 3: return False
                    h_len = 0
                # Vertical check
                if self.grid[j][i] == 0:
                    v_len += 1
                else:
                    if 0 < v_len < 3: return False
                    v_len = 0
            if 0 < h_len < 3 or 0 < v_len < 3:
                return False
        return True

if __name__ == '__main__':
    grid_size = 8
    counter = CrosswordGridCounter(grid_size)
    num_grids = counter.count_grids()
    # The prompt asked to output each number in the final equation.
    # As there is no equation, we will just print the final result clearly.
    print(f"For an {grid_size}x{grid_size} puzzle with the specified rules:")
    print(f"The number of possible grids is: {num_grids}")