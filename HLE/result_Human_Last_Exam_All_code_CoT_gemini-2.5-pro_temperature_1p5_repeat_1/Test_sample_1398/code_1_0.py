import collections

def solve():
    """
    Solves the grid transformation puzzle.
    The rule is that a cell becomes 1 if it has 2 or 3 neighbors, and 0 otherwise.
    The original value of the cell does not matter.
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

    # Iterate through each cell of the grid
    for r in range(rows):
        for c in range(cols):
            # Count the number of '1's in the 8 neighboring cells
            neighbor_sum = 0
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    # Skip the cell itself
                    if dr == 0 and dc == 0:
                        continue
                    
                    # Get neighbor coordinates
                    nr, nc = r + dr, c + dc
                    
                    # Check if the neighbor is within the grid boundaries
                    if 0 <= nr < rows and 0 <= nc < cols:
                        neighbor_sum += input_grid[nr][nc]
            
            # Apply the rule: if neighbor count is 2 or 3, the cell is 1
            if neighbor_sum == 2 or neighbor_sum == 3:
                output_grid[r][c] = 1
            else:
                output_grid[r][c] = 0

    # Print the final output grid
    final_string = ""
    for r in range(rows):
        row_str = " ".join(map(str, output_grid[r]))
        print(f"[{row_str}]")
        final_string += "".join(map(str, output_grid[r]))

    print("\nFinal answer string:")
    print(final_string)

solve()