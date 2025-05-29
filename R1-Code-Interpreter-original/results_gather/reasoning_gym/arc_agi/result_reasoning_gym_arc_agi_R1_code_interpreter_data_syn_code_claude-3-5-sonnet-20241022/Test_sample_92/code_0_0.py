def expand_pattern(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output = [[3 for _ in range(cols)] for _ in range(rows)]
    
    # Helper function to check if a position should be filled
    def should_fill(r, c, num, positions):
        count = 0
        for pos in positions:
            pr, pc = pos
            if 0 <= pr < rows and 0 <= pc < cols and input_grid[pr][pc] == num:
                count += 1
        return count >= 1

    # Process each position in the grid
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] != 3:
                val = input_grid[i][j]
                
                # Define expansion pattern based on the number
                if val in [5, 6, 7, 8, 9]:
                    # These numbers tend to expand in a cross or diamond pattern
                    positions = [(i,j)]
                    for di, dj in [(0,1), (1,0), (-1,0), (0,-1)]:
                        ni, nj = i + di, j + dj
                        if 0 <= ni < rows and 0 <= nj < cols:
                            positions.append((ni, nj))
                    
                    # Fill the pattern
                    for pos in positions:
                        pr, pc = pos
                        if should_fill(pr, pc, val, [(i,j)]):
                            output[pr][pc] = val
                
                elif val in [1]:
                    # 1 tends to expand in small clusters
                    positions = [(i,j)]
                    for di, dj in [(0,1), (1,0), (-1,0), (0,-1)]:
                        ni, nj = i + di, j + dj
                        if 0 <= ni < rows and 0 <= nj < cols:
                            positions.append((ni, nj))
                    
                    for pos in positions:
                        pr, pc = pos
                        if should_fill(pr, pc, val, [(i,j)]):
                            output[pr][pc] = val

    # Copy original numbers
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] != 3:
                output[i][j] = input_grid[i][j]
                
    return output

# Parse input grid
input_str = """3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 7 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 7 7 7 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 9 7 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 1 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 1 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 6 1 1 1 3 3 3 3 3 3 5 9 3 3 3 3 3 3 3
3 3 3 3 3 7 5 3 3 3 3 3 3 3 3 7 6 6 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 6 3 6 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 6 6 6 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 6 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3"""

input_grid = [[int(x) for x in row.split()] for row in input_str.split('\n')]
output_grid = expand_pattern(input_grid)

# Format output
output_str = '\n'.join(' '.join(str(x) for x in row) for row in output_grid)
print(f"<<<{output_str}>>>")