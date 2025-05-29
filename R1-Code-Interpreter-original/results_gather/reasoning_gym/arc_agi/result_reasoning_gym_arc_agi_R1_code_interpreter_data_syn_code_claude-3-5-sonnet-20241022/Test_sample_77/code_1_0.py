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
        for i in range(2, 8):
            if i + 3 < 8:  # Check if we have room for a sequence of 4
                # Check if we have a sequence of identical numbers
                value = grid[i][j]
                seq_len = 1
                while i + seq_len < 8 and grid[i+seq_len][j] == value:
                    seq_len += 1
                
                # If sequence is at least 4 numbers long, replace with bottom half numbers
                if seq_len >= 4:
                    for k in range(seq_len):
                        output[i+k][j] = grid[13-k][j]
        
        # Check for sequences in bottom half (rows 8-13)
        for i in range(8, 14):
            if i + 3 < 14:  # Check if we have room for a sequence of 4
                # Check if we have a sequence of identical numbers
                value = grid[i][j]
                seq_len = 1
                while i + seq_len < 14 and grid[i+seq_len][j] == value:
                    seq_len += 1
                
                # If sequence is at least 4 numbers long, replace with top half numbers
                if seq_len >= 4:
                    for k in range(seq_len):
                        output[i+k][j] = grid[5-k][j]
    
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