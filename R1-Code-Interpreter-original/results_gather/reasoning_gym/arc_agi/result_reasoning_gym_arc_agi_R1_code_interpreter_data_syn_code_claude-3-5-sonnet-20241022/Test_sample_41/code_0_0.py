def process_grid(grid):
    n = len(grid)
    output = [row[:] for row in grid]
    
    # Process rows for sequences of 3's
    for i in range(n):
        for j in range(n):
            # Check for horizontal sequences of 3's
            if j < n-3 and all(grid[i][j+k] == 3 for k in range(4)):
                # Get value from above or below
                if i > 0:
                    val = grid[i-1][j]
                    for k in range(4):
                        output[i][j+k] = val
                elif i < n-1:
                    val = grid[i+1][j]
                    for k in range(4):
                        output[i][j+k] = val
            
            # Check for vertical sequences of 3's
            if i < n-3 and all(grid[i+k][j] == 3 for k in range(4)):
                # Get value from left or right
                if j > 0:
                    val = grid[i][j-1]
                    for k in range(4):
                        output[i+k][j] = val
                elif j < n-1:
                    val = grid[i][j+1]
                    for k in range(4):
                        output[i+k][j] = val

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