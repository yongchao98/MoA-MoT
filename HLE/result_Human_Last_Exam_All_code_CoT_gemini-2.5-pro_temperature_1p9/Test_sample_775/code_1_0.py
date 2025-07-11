import numpy as np
from collections import deque

def solve_puzzle():
    """
    Analyzes the intersection of two specially constructed sets A and B.

    The sets are defined on a grid and designed to test the connectivity
    of their intersection. This code provides a computational example for a
    topology problem.
    """
    # 1. Define the grid
    grid_size = 100
    grid_A = np.zeros((grid_size, grid_size), dtype=bool)
    grid_B = np.zeros((grid_size, grid_size), dtype=bool)

    # 2. Define Set A
    # A consists of a top block, a bottom block, and a thin vertical spine
    # on the left to connect them.
    bottom_block_end = 33
    top_block_start = 66
    spine_width = 2

    # Bottom block
    grid_A[0:bottom_block_end, :] = True
    # Top block
    grid_A[top_block_start:, :] = True
    # Left spine (to make A connected)
    grid_A[:, 0:spine_width] = True

    # 3. Define Set B
    # B is the middle block. Since grid_A and grid_B are 'closed' sets on the
    # grid, their boundaries will overlap.
    grid_B[bottom_block_end-1:top_block_start+1, :] = True
    
    # 4. Compute the intersection
    # A cell is in the intersection if it is True in both grids.
    grid_intersection = np.logical_and(grid_A, grid_B)

    # 5. Count connected components in the intersection
    def count_components(grid):
        """
        Counts connected components in a boolean grid using BFS.
        A component is a group of 'True' cells connected 8-directionally.
        """
        if not np.any(grid):
            return 0
        
        rows, cols = grid.shape
        visited = np.zeros((rows, cols), dtype=bool)
        num_components = 0
        
        for r in range(rows):
            for c in range(cols):
                if grid[r, c] and not visited[r, c]:
                    num_components += 1
                    queue = deque([(r, c)])
                    visited[r, c] = True
                    
                    while queue:
                        row, col = queue.popleft()
                        # Check 8 neighbors
                        for dr in [-1, 0, 1]:
                            for dc in [-1, 0, 1]:
                                if dr == 0 and dc == 0:
                                    continue
                                
                                nr, nc = row + dr, col + dc
                                
                                if 0 <= nr < rows and 0 <= nc < cols and \
                                   grid[nr, nc] and not visited[nr, nc]:
                                    visited[nr, nc] = True
                                    queue.append((nr, nc))
        return num_components

    num_A_components = count_components(grid_A)
    num_B_components = count_components(grid_B)
    num_union_components = count_components(np.logical_or(grid_A, grid_B))
    num_intersection_components = count_components(grid_intersection)
    
    print("This code computationally investigates a tricky case.")
    print("The properties of the constructed sets on the grid are:")
    print(f"Number of components in A: {num_A_components}")
    print(f"Number of components in B: {num_B_components}")
    print(f"Number of components in A U B: {num_union_components}")
    print("\n--- Final Equation ---")
    print(f"Largest number of components of the intersection = {num_intersection_components}")


solve_puzzle()