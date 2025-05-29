def analyze_pattern(input_grid):
    # Convert input grid to 2D array
    grid = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    rows = len(grid)
    cols = len(grid[0])
    
    # Find all horizontal lines
    line_rows = []
    for i, row in enumerate(grid):
        if len(set(row)) == 1:  # All numbers in row are the same
            line_rows.append(i)
    
    # Find the most common number in non-line rows
    non_line_nums = []
    for i, row in enumerate(grid):
        if i not in line_rows:
            non_line_nums.extend(row)
    
    from collections import Counter
    most_common = Counter(non_line_nums).most_common(1)[0][0]
    
    # Calculate distances between lines
    distances = []
    prev = -1
    for line in line_rows:
        if prev >= 0:
            distances.append(line - prev - 1)
        prev = line
    
    print(f"Input dimensions: {rows}x{cols}")
    print(f"Line rows: {line_rows}")
    print(f"Most common non-line number: {most_common}")
    print(f"Number of horizontal lines: {len(line_rows)}")
    print(f"Distances between lines: {distances}")

# Test with the full input
test_input = """2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2"""

analyze_pattern(test_input)