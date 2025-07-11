# The size of the grid
n = 4
m = 4

# Total number of cells in the grid
num_cells = n * m

# Total number of possible assignments (2^16 for a 4x4 grid)
total_assignments = 1 << num_cells

# Counter for valid assignments
valid_count = 0

# Iterate through each possible assignment, represented by an integer from 0 to 2^16 - 1
for i in range(total_assignments):
    # Create the 4x4 grid for the current assignment 'i'
    # We use bitwise operations to map the integer 'i' to the grid cells.
    # The bit at position (r * m + c) corresponds to the cell (r, c).
    grid = [[(i >> (r * m + c)) & 1 for c in range(m)] for r in range(n)]
    
    # Assume the grid is valid until a violation is found
    is_valid = True
    
    # Check the validity of the current grid
    for r in range(n):
        for c in range(m):
            # If the current cell is 1, check its right and bottom neighbors
            if grid[r][c] == 1:
                # Check for horizontal violation (with the cell to the right)
                if c + 1 < m and grid[r][c+1] == 1:
                    is_valid = False
                    break
                # Check for vertical violation (with the cell below)
                if r + 1 < n and grid[r+1][c] == 1:
                    is_valid = False
                    break
        if not is_valid:
            break
            
    # If the grid is valid after all checks, increment the counter
    if is_valid:
        valid_count += 1

# Output the results of the calculation
print(f"Grid dimensions: n = {n}, m = {m}")
print(f"Total number of cells: {n} * {m} = {num_cells}")
print(f"Total number of assignments to check: 2^{num_cells} = {total_assignments}")
print(f"Final count of valid assignments: {valid_count}")