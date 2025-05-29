def process_grid(grid):
    n = len(grid)
    output = [row[:] for row in grid]
    
    # Process middle section of the grid
    for i in range(n):
        for j in range(n):
            # If we find a 3, look for pattern of 3's
            if grid[i][j] == 3:
                # Check if it's part of a horizontal sequence
                if j+3 < n and all(grid[i][k] == 3 for k in range(j, j+4)):
                    # Look for the value to replace with from row above
                    if i > 0:
                        for k in range(j, j+4):
                            output[i][k] = grid[i-1][k]
                
                # Check if it's part of a vertical sequence
                if i+3 < n and all(grid[k][j] == 3 for k in range(i, i+4)):
                    # Look for the value to replace with from column to the left
                    if j > 0:
                        for k in range(i, i+4):
                            output[k][j] = grid[k][j-1]

    return output

# Parse input grid
input_grid = [
    [4,4,7,5,2,7,9,8,8,9,7,2,5,7,4,4],
    [4,5,5,7,7,7,8,8,8,8,7,7,7,5,5,4],
    [7,5,8,2,9,8,8,3,3,3,8,9,2,8,5,7],
    [5,7,2,5,8,8,5,3,3,3,8,8,5,2,7,5],
    [2,7,3,3,3,3,8,3,3,3,4,7,8,9,7,2],
    [7,7,3,3,3,3,8,0,0,8,4,4,8,8,7,7],
    [9,8,3,3,3,3,0,0,0,0,8,8,5,8,8,9],
    [8,8,3,3,3,3,0,4,4,0,0,8,5,5,8,8],
    [8,8,5,5,8,0,0,4,4,0,0,8,5,5,8,8],
    [9,8,8,5,8,8,0,0,0,0,8,8,5,8,8,9],
    [7,7,8,8,4,4,8,0,0,8,4,4,8,8,7,7],
    [2,7,9,8,7,4,8,8,8,8,4,7,8,9,7,2],
    [5,7,2,5,8,8,5,5,5,5,8,8,5,2,7,5],
    [7,5,8,2,9,8,8,5,5,8,8,9,2,8,5,7],
    [4,5,5,7,7,7,8,8,8,8,7,7,7,5,5,4],
    [4,4,7,5,2,7,9,8,8,9,7,2,5,7,4,4]
]

result = process_grid(input_grid)

# Print the result in the required format
print("<<<")
for row in result:
    print(" ".join(map(str, row)))
print(">>>")