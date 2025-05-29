def analyze_grid(grid):
    # Convert string grid to 2D array of numbers
    rows = len(grid)
    cols = len(grid[0].split())
    num_grid = [[int(x) for x in row.split()] for row in grid]
    
    # Find the central region
    start_row = rows // 3
    end_row = rows - rows // 3
    start_col = cols // 3
    end_col = cols - cols // 3
    
    # Extract central region numbers
    central_numbers = []
    for i in range(start_row, end_row):
        for j in range(start_col, end_col):
            central_numbers.append(num_grid[i][j])
    
    # Count frequency of numbers in central region
    from collections import Counter
    number_freq = Counter(central_numbers)
    most_common = number_freq.most_common(2)
    
    print("Central region analysis:")
    print(f"Dimensions: {end_row - start_row} x {end_col - start_col}")
    print("Most common numbers:", most_common)
    print("\nCentral region pattern:")
    for i in range(start_row, min(end_row, start_row + 5)):
        central_row = num_grid[i][start_col:end_col]
        print(f"Row {i}:", central_row)

# Test input
test_input = """9 4 3 7 0 2 3 3 5 7 4 1 7 8 5 7 9 8 6
8 4 8 5 9 9 4 7 4 8 0 8 7 1 9 6 6 6 9
9 3 5 1 5 7 4 9 0 8 1 8 4 2 9 5 9 5 2
0 8 3 3 3 3 3 3 3 3 3 3 8 6 2 7 4 4 8
9 0 3 3 3 3 3 3 3 3 3 3 7 4 4 1 5 9 3
7 8 3 3 3 3 3 6 3 3 3 3 3 6 3 6 9 3 8
1 5 3 3 3 3 3 3 3 3 3 3 6 7 1 7 6 4 6
0 6 3 3 3 3 3 3 3 3 3 3 3 4 5 7 2 1 7
1 4 3 3 3 3 3 3 3 3 3 3 6 5 5 3 6 7 2
5 5 3 3 6 3 3 3 3 3 3 3 5 7 5 7 8 9 3
6 2 3 3 3 3 3 3 3 3 3 3 9 8 2 0 0 3 0
1 9 8 8 5 3 3 7 7 4 4 4 7 8 5 6 8 0 5""".split('\n')

analyze_grid(test_input)