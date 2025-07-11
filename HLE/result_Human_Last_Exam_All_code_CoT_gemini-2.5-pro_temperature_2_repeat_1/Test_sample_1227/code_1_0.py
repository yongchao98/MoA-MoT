import sys

# The recursive search can be deep, so we increase the recursion limit.
sys.setrecursionlimit(2000)

class CrosswordGridCounter:
    """
    A class to find and count valid crossword grids based on a set of rules.
    """
    def __init__(self, size):
        if size % 2 != 0:
            raise ValueError("Size must be an even number for this implementation.")
        self.N = size
        self.grid = [[-1] * self.N for _ in range(self.N)]
        self.count = 0

    def solve(self):
        """
        Starts the backtracking search to find all valid grid configurations.
        """
        # We need to fill N*N/2 cells for a grid of size N.
        # For 8x8, this is 32 cells (the top 4 rows).
        self.backtrack(0)
        return self.count

    def check_connectivity(self):
        """
        Checks if all white squares (1) are part of a single connected component.
        Returns False if there are no white squares or they are not connected.
        """
        white_squares = []
        for r in range(self.N):
            for c in range(self.N):
                if self.grid[r][c] == 1:
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
                if 0 <= nr < self.N and 0 <= nc < self.N and self.grid[nr][nc] == 1 and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))
                    
        return len(visited) == len(white_squares)

    def check_word_lengths(self):
        """
        Checks that there are no words of length 1 or 2.
        A "word" is a contiguous horizontal or vertical run of white squares.
        """
        for r in range(self.N): # Check rows
            length = 0
            for c in range(self.N):
                if self.grid[r][c] == 1:
                    length += 1
                else:
                    if 0 < length < 3: return False
                    length = 0
            if 0 < length < 3: return False
            
        for c in range(self.N): # Check columns
            length = 0
            for r in range(self.N):
                if self.grid[r][c] == 1:
                    length += 1
                else:
                    if 0 < length < 3: return False
                    length = 0
            if 0 < length < 3: return False
            
        return True

    def check_cheaters_pruning(self, r_placed, c_placed):
        """
        Checks for violations of the "no cheater" rule for pruning purposes.
        The rule is: no 2x2 block can have 3 or more black squares.
        This check is limited to 2x2 blocks that have been fully determined.
        """
        points_to_check = {(r_placed, c_placed), (self.N - 1 - r_placed, self.N - 1 - c_placed)}

        for r_center, c_center in points_to_check:
            # Iterate through the 4 possible 2x2 blocks containing (r_center, c_center)
            # (dr_start, dc_start) is the offset of the top-left corner from the center point
            for dr_start in [-1, 0]:
                for dc_start in [-1, 0]:
                    r_start, c_start = r_center + dr_start, c_center + dc_start

                    if not (0 <= r_start < self.N - 1 and 0 <= c_start < self.N - 1):
                        continue
                    
                    black_count = 0
                    all_determined = True
                    for i in range(r_start, r_start + 2):
                        for j in range(c_start, c_start + 2):
                            val = self.grid[i][j]
                            if val == -1:
                                all_determined = False; break
                            elif val == 0:
                                black_count += 1
                        if not all_determined: break
                    
                    if all_determined and black_count >= 3:
                        return False
        return True

    def backtrack(self, k):
        """
        The main recursive backtracking function.
        k: the index of the cell to be filled (0 to 31 for an 8x8 grid).
        """
        # Base case: The top half of the grid is filled. The full grid is now determined.
        if k == (self.N * self.N) // 2:
            if self.check_word_lengths() and self.check_connectivity():
                self.count += 1
            return

        r, c = k // self.N, k % self.N
        
        # Choice 1: Place a white square.
        self.grid[r][c] = 1
        self.grid[self.N - 1 - r][self.N - 1 - c] = 1
        self.backtrack(k + 1)
        
        # Choice 2: Place a black square.
        self.grid[r][c] = 0
        self.grid[self.N - 1 - r][self.N - 1 - c] = 0
        # Prune this branch if the placement of black squares creates an invalid state.
        if self.check_cheaters_pruning(r, c):
            self.backtrack(k + 1)

        # Backtrack: reset the cell to undecided for other search paths.
        self.grid[r][c] = -1
        self.grid[self.N - 1 - r][self.N - 1 - c] = -1

if __name__ == '__main__':
    # Note: This is a computationally intensive task and may take a very long time to complete.
    grid_size = 8
    counter = CrosswordGridCounter(grid_size)
    num_grids = counter.solve()
    print(f"For a {grid_size}x{grid_size} puzzle with the specified rules, the number of possible grids is:")
    print(num_grids)