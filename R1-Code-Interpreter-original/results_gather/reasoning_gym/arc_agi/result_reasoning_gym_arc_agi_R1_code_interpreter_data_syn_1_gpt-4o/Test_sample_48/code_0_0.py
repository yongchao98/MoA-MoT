def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid

    def flood_fill(r, c, visited):
        # Perform a flood fill to find all connected '1's
        if r < 0 or r >= rows or c < 0 or c >= cols:
            return []
        if (r, c) in visited or input_grid[r][c] != 1:
            return []
        
        visited.add((r, c))
        block = [(r, c)]
        
        # Check all 4 directions
        block += flood_fill(r + 1, c, visited)
        block += flood_fill(r - 1, c, visited)
        block += flood_fill(r, c + 1, visited)
        block += flood_fill(r, c - 1, visited)
        
        return block

    visited = set()
    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] == 1 and (r, c) not in visited:
                block = flood_fill(r, c, visited)
                if len(block) > 2:  # Only replace if the block is larger than 2
                    for br, bc in block:
                        output_grid[br][bc] = 5

    return output_grid

# Test input grid
input_grid = [
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 1],
    [6, 6, 6, 1, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 1, 6, 6, 6, 1, 1, 1, 1, 1, 6, 6, 6],
    [6, 6, 6, 1, 6, 6, 6, 1, 6, 6, 1, 6, 6, 6, 6],
    [6, 6, 6, 1, 6, 6, 6, 1, 6, 6, 1, 6, 6, 6, 6],
    [1, 6, 6, 1, 6, 6, 6, 1, 1, 1, 1, 6, 6, 6, 6],
    [1, 1, 1, 1, 6, 6, 6, 1, 6, 6, 6, 6, 6, 6, 6],
    [1, 6, 6, 1, 6, 6, 6, 1, 6, 6, 6, 1, 1, 1, 1],
    [1, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 6, 6, 1],
    [1, 1, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 6, 6, 1],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 1],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 1, 1, 6, 6, 6, 1, 1, 1, 1, 1, 1, 6, 6, 6],
    [6, 1, 6, 6, 6, 6, 1, 6, 6, 6, 6, 1, 6, 6, 6],
    [6, 1, 6, 6, 6, 6, 1, 6, 6, 6, 6, 1, 6, 6, 6],
    [6, 1, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 6, 6, 6]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))