def convert_to_grid(input_str):
    rows = []
    numbers = input_str.strip().split()
    for i in range(0, len(numbers), 16):
        rows.append([int(x) for x in numbers[i:i+16]])
    return rows

def process_grid(grid):
    n = 16
    output = [row[:] for row in grid]
    
    # Check vertical sequences in rows 2-5 and 10-13
    for j in range(n):
        # Check upper part (rows 2-5)
        if (len(set([grid[2][j], grid[3][j], grid[4][j], grid[5][j]])) == 1):
            # Copy from corresponding lower part (rows 10-13)
            output[2][j] = grid[13][j]
            output[3][j] = grid[12][j]
            output[4][j] = grid[11][j]
            output[5][j] = grid[10][j]
        
        # Check lower part (rows 10-13)
        if (len(set([grid[10][j], grid[11][j], grid[12][j], grid[13][j]])) == 1):
            # Copy from corresponding upper part (rows 2-5)
            output[13][j] = grid[2][j]
            output[12][j] = grid[3][j]
            output[11][j] = grid[4][j]
            output[10][j] = grid[5][j]
    
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