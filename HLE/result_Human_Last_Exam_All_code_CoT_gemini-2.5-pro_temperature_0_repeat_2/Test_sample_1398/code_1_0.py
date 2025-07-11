import collections

def solve():
    """
    This function solves the grid transformation puzzle by applying a specific cellular automaton rule.
    The rule is: a cell in the output grid is 1 if its corresponding cell in the input grid
    has exactly 2 or 3 live neighbors (value 1) in its 8-cell Moore neighborhood. Otherwise, it's 0.
    """
    input_grid = [
        [0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 0, 1],
        [1, 0, 1, 1, 0]
    ]

    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [[0 for _ in range(cols)] for _ in range(rows)]

    for r in range(rows):
        for c in range(cols):
            live_neighbors = 0
            # Iterate over the 8 neighbors (Moore neighborhood)
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    if dr == 0 and dc == 0:
                        continue
                    
                    nr, nc = r + dr, c + dc

                    # Check if the neighbor is within the grid boundaries
                    if 0 <= nr < rows and 0 <= nc < cols:
                        live_neighbors += input_grid[nr][nc]
            
            # Apply the rule: cell is 1 if it has 2 or 3 live neighbors
            if live_neighbors == 2 or live_neighbors == 3:
                output_grid[r][c] = 1
            else:
                output_grid[r][c] = 0

    # Print the resulting grid
    for row in output_grid:
        print(' '.join(map(str, row)))

solve()