def find_valid_structure(grid, r, c):
    rows, cols = len(grid), len(grid[0])
    
    # Check if position is center of a 3x3 pattern with 9 in center and 6s around
    def is_valid_center(row, col):
        if row < 1 or row >= rows-1 or col < 1 or col >= cols-1:
            return False
        if grid[row][col] != 9:
            return False
        # Check surrounding positions should be 6s
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                if grid[row+dr][col+dc] != 6:
                    return False
        return True
    
    # Find all connected positions that are part of valid structures
    valid_positions = set()
    to_check = [(r, c)]
    checked = set()
    
    while to_check:
        curr_r, curr_c = to_check.pop(0)
        if (curr_r, curr_c) in checked:
            continue
            
        checked.add((curr_r, curr_c))
        
        # If current position is center of valid pattern
        if is_valid_center(curr_r, curr_c):
            # Add all positions in the 3x3 grid
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    valid_positions.add((curr_r+dr, curr_c+dc))
                    # Add adjacent positions to check
                    if grid[curr_r+dr][curr_c+dc] in [6, 9]:
                        to_check.append((curr_r+dr, curr_c+dc))
                        
        # If current position is a 6 adjacent to a valid pattern
        elif grid[curr_r][curr_c] == 6:
            # Check neighbors for valid centers
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    new_r, new_c = curr_r+dr, curr_c+dc
                    if (0 <= new_r < rows and 0 <= new_c < cols and 
                        grid[new_r][new_c] == 9):
                        to_check.append((new_r, new_c))
    
    return valid_positions

def transform_grid(input_grid):
    rows, cols = len(input_grid), len(input_grid[0])
    output_grid = [[8 for _ in range(cols)] for _ in range(rows)]
    
    # Find all valid structures connected to bottom row
    valid_positions = set()
    for c in range(cols):
        if input_grid[rows-1][c] in [6, 9]:
            valid_positions.update(find_valid_structure(input_grid, rows-1, c))
    
    # Build output grid
    for r in range(rows):
        for c in range(cols):
            if (r, c) in valid_positions:
                output_grid[r][c] = input_grid[r][c]
    
    return output_grid

# Test input
test_input = [
    [8,8,8,8,8,8,8,8,6,6,6,8],
    [8,6,6,6,6,6,6,8,6,6,6,8],
    [8,6,9,6,6,6,6,8,8,8,8,8],
    [8,6,6,6,6,9,6,8,6,6,6,6],
    [8,6,6,6,6,6,6,8,6,6,6,6],
    [8,8,8,8,8,8,8,8,6,6,9,6],
    [8,8,8,8,6,6,6,8,6,6,6,6],
    [8,8,8,8,6,6,9,8,8,8,8,8],
    [8,8,8,8,6,6,6,8,8,8,8,8],
    [6,6,6,8,9,6,6,8,8,8,8,8],
    [9,6,6,8,8,8,8,8,6,9,6,9],
    [6,6,6,8,8,8,8,8,6,6,6,6],
    [8,8,8,8,8,8,8,8,6,9,6,6]
]

output = transform_grid(test_input)
print("Output grid:")
for row in output:
    print(' '.join(map(str, row)))