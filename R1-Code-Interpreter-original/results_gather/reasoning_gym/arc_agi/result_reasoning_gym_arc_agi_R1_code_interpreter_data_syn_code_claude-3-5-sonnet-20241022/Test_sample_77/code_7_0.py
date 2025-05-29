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
        # Check for sequences in top half (rows 2-7)
        for i in range(2, 7):
            if (i + 3 < 8 and  # Ensure we have room for a sequence of 4
                len(set([grid[i][j], grid[i+1][j], grid[i+2][j], grid[i+3][j]])) == 1):
                # Found a sequence, replace with numbers from bottom half
                mirror_start = 13
                for k in range(i, i+4):
                    output[k][j] = grid[mirror_start-k+2][j]
        
        # Check for sequences in bottom half (rows 8-13)
        for i in range(8, 14):
            if (i + 3 < 14 and  # Ensure we have room for a sequence of 4
                len(set([grid[i][j], grid[i+1][j], grid[i+2][j], grid[i+3][j]])) == 1):
                # Found a sequence, replace with numbers from top half
                mirror_start = 5
                for k in range(i, i+4):
                    output[k][j] = grid[mirror_start-(k-8)][j]
    
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