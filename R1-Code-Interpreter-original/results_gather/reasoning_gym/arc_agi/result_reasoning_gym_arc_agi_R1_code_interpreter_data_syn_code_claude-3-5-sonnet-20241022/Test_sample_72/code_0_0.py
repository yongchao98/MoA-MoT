def analyze_grid(input_grid):
    # Convert string grid to 2D array of numbers
    grid = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    
    # Find the horizontal lines (rows where all numbers are the same)
    line_rows = []
    line_numbers = set()
    for i, row in enumerate(grid):
        if len(set(row)) == 1:  # All numbers in row are the same
            line_rows.append(i)
            line_numbers.add(row[0])
    
    # Find the most common number that's not in a line
    numbers = []
    for i, row in enumerate(grid):
        if i not in line_rows:
            numbers.extend(row)
    
    if numbers:
        from collections import Counter
        most_common = Counter(numbers).most_common(1)[0][0]
    else:
        most_common = None
    
    # Count sections between lines
    sections = len(line_rows) + 1
    
    print(f"Line numbers: {line_numbers}")
    print(f"Most common non-line number: {most_common}")
    print(f"Number of sections: {sections}")
    
    # For the test input
    test_input = """2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2"""
    
    # Get dimensions
    rows = len(grid)
    cols = len(grid[0])
    print(f"Grid dimensions: {rows}x{cols}")

# Test with the first few rows of the test input
analyze_grid("""2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
2 7 2 2 2 2 7 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2""")