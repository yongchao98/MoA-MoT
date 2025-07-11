def solve():
    """
    This function applies a cellular automaton rule (B23/S23) to a given grid.
    - A cell is BORN (0->1) if it has 2 or 3 live neighbors.
    - A cell SURVIVES (1->1) if it has 2 or 3 live neighbors.
    - Otherwise, a cell DIES (1->0 or 0->0).
    The neighborhood is the 8 surrounding cells (Moore neighborhood), without wrap-around.
    """
    
    # The sixth test input grid
    grid = [
        [0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 0, 1],
        [1, 0, 1, 1, 0]
    ]

    rows = len(grid)
    cols = len(grid[0])
    output_grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Iterate through each cell of the grid to calculate the new state
    for r in range(rows):
        for c in range(cols):
            # Count live Moore neighbors (8 directions)
            live_neighbors = 0
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    # Skip the cell itself
                    if dr == 0 and dc == 0:
                        continue
                    
                    nr, nc = r + dr, c + dc
                    
                    # Check if the neighbor is within grid boundaries
                    if 0 <= nr < rows and 0 <= nc < cols:
                        if grid[nr][nc] == 1:
                            live_neighbors += 1

            # Apply the B23/S23 rule
            current_value = grid[r][c]
            
            if current_value == 1:  # Live cell
                if live_neighbors == 2 or live_neighbors == 3:
                    output_grid[r][c] = 1  # Survival
                else:
                    output_grid[r][c] = 0  # Death by underpopulation or overpopulation
            else:  # Dead cell (current_value == 0)
                if live_neighbors == 2 or live_neighbors == 3:
                    output_grid[r][c] = 1  # Birth
                else:
                    output_grid[r][c] = 0  # Stays dead

    # Print the final output grid row by row
    for row in output_grid:
        print(' '.join(map(str, row)))

solve()