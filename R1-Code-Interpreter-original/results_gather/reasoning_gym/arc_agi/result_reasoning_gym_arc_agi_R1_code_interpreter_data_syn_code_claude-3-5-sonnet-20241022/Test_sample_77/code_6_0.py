def convert_to_grid(input_str):
    rows = []
    numbers = input_str.strip().split()
    for i in range(0, len(numbers), 16):
        rows.append([int(x) for x in numbers[i:i+16]])
    return rows

def process_grid(grid):
    n = 16
    output = [row[:] for row in grid]
    
    # Process each column
    for j in range(n):
        # Check for sequences in both top and bottom halves simultaneously
        for i in range(2, 6):  # Check rows 2-5 (top half)
            bottom_i = n - i - 1  # Corresponding row in bottom half
            
            # Check if we have a sequence of 4 identical numbers in either half
            if (grid[i][j] == grid[i+1][j] == grid[i+2][j] == grid[i+3][j]):
                # Replace top sequence with bottom numbers
                output[i][j] = grid[bottom_i][j]
                output[i+1][j] = grid[bottom_i-1][j]
                output[i+2][j] = grid[bottom_i-2][j]
                output[i+3][j] = grid[bottom_i-3][j]
            
            if (grid[bottom_i][j] == grid[bottom_i-1][j] == grid[bottom_i-2][j] == grid[bottom_i-3][j]):
                # Replace bottom sequence with top numbers
                output[bottom_i][j] = grid[i][j]
                output[bottom_i-1][j] = grid[i+1][j]
                output[bottom_i-2][j] = grid[i+2][j]
                output[bottom_i-3][j] = grid[i+3][j]
    
    return output

input_str = """3 3 1 8 0 1 5 4 4 5 1 0 8 1 3 3
3 8 8 1 1 1 4 4 4 4 1 1 1 8 8 3
1 8 4 0 5 4 4 2 2 2 4 5 0 4 8 1
8 1 0 8 4 4 8 2 2 2 4 4 8 0 1 8
0 1 2 2 2 2 4 2 2 2 3 1 4 5 1 0
1 1 2 2 2 2 4 7 7 4 3 3 4 4 1 1
5 4 2 2 2 2 7 7 7 7 4 4 8 4 4 5
4 4 2 2 2 2 7 3 3 7 7 4 8 8 4 4
4 4 8 8 4 7 7 3 3 7 7 4 8 8 4 4
5 4 4 8 4 4 7 7 7 7 4 4 8 4 4 5
1 1 4 4 3 3 4 7 7 4 3 3 4 4 1 1
0 1 5 4 1 3 4 4 4 4 3 1 4 5 1 0
8 1 0 8 4 4 8 8 8 8 4 4 8 0 1 8
1 8 4 0 5 4 4 8 8 4 4 5 0 4 8 1
3 8 8 1 1 1 4 4 4 4 1 1 1 8 8 3
3 3 1 8 0 1 5 4 4 5 1 0 8 1 3 3"""

grid = convert_to_grid(input_str)
result = process_grid(grid)

for row in result:
    print(' '.join(map(str, row)))