def count_rectangles(grid):
    if not grid or not grid[0]:
        return 0
    
    rows, cols = len(grid), len(grid[0])
    count = 0
    
    # Store horizontal lines: (start_col, end_col) for each row
    h_lines = []
    for i in range(rows):
        lines = []
        start = -1
        for j in range(cols):
            if grid[i][j] in ['#', '█']:
                if start == -1:
                    start = j
            elif start != -1:
                lines.append((start, j-1))
                start = -1
        if start != -1:
            lines.append((start, cols-1))
        h_lines.append(lines)
    
    # For each row, look for matching horizontal lines below
    for i in range(rows):
        for line1 in h_lines[i]:
            start_col, end_col = line1
            # Look for matching bottom line
            for k in range(i+1, rows):
                for line2 in h_lines[k]:
                    if line2 == line1:  # Same horizontal span
                        # Verify vertical lines
                        valid = True
                        for col in [start_col, end_col]:
                            # Check if vertical lines exist
                            for row in range(i, k+1):
                                if grid[row][col] not in ['#', '█']:
                                    valid = False
                                    break
                            if not valid:
                                break
                        if valid:
                            count += 1
    
    return count

# Read the grid
grid = []
while True:
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

print(count_rectangles(grid))