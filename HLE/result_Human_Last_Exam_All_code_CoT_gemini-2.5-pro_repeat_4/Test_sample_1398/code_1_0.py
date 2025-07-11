import numpy as np

def solve():
    """
    Applies a discovered cellular automaton rule to a given grid.
    The rule is B2567/S23 on a Moore neighborhood.
    - A cell is "born" (0->1) if it has 2, 5, 6, or 7 neighbors.
    - A cell "survives" (1->1) if it has 2 or 3 neighbors.
    """
    
    # The 6th input grid provided in the problem
    input_grid_data = [
        [0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 0, 1],
        [1, 0, 1, 1, 0]
    ]
    
    input_grid = np.array(input_grid_data)
    rows, cols = input_grid.shape
    output_grid = np.zeros_like(input_grid)

    for r in range(rows):
        for c in range(cols):
            # Calculate the sum of the 8 Moore neighbors
            # We treat the grid as being padded with 0s at the boundaries
            s_moore = 0
            for i in range(max(0, r - 1), min(rows, r + 2)):
                for j in range(max(0, c - 1), min(cols, c + 2)):
                    if i == r and j == c:
                        continue
                    s_moore += input_grid[i, j]
            
            # Apply the B2567/S23 rule
            current_value = input_grid[r, c]
            new_value = 0
            if current_value == 0:  # Birth rule
                if s_moore in [2, 5, 6, 7]:
                    new_value = 1
            else:  # Survival rule (current_value == 1)
                if s_moore in [2, 3]:
                    new_value = 1
            
            output_grid[r, c] = new_value

    # Print the final output grid, row by row
    for row in output_grid:
        # The prompt asks to "output each number in the final equation"
        # We interpret this as printing the final state of each number/cell.
        print(' '.join(map(str, row)))

solve()