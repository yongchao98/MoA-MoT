def find_central_region(grid):
    # Convert string grid to 2D array of numbers
    rows = len(grid)
    cols = len(grid[0].split())
    num_grid = [[int(x) for x in row.split()] for row in grid]
    
    # Find the most common number in the central region
    central_numbers = []
    for i in range(rows):
        for j in range(cols):
            if len(set(row.split()) for row in grid) > 1:  # if not all same numbers
                central_numbers.append(num_grid[i][j])
    
    # Print analysis of the central region
    print("Central region analysis:")
    print("Most common numbers:", sorted(list(set(central_numbers))))
    
    # Analyze the first example
    print("\nAnalyzing pattern in central region:")
    for i in range(4, 13):  # Focus on middle rows
        middle_nums = num_grid[i][6:16]  # Focus on middle columns
        print(f"Row {i}:", middle_nums)

# Example 1 input
example1 = """1 2 7 3 0 3 2 0 8 5 2 4 8 9 7 1 6 4 9 5 1 6 9 8 4 9 9
6 7 3 2 1 4 7 0 1 1 0 7 4 2 4 3 0 7 0 4 3 5 1 7 6 2 2
3 3 7 8 3 8 0 6 1 0 3 4 3 9 2 7 6 4 7 1 9 4 5 7 8 4 2
9 1 0 4 9 5 3 1 6 7 6 1 9 4 7 6 4 7 7 3 9 2 1 4 7 4 0
2 5 9 4 7 2 4 4 4 4 4 4 4 4 4 4 2 9 0 7 1 3 5 6 7 0 0
6 3 0 4 4 9 4 4 4 4 4 4 4 4 4 4 0 8 1 9 7 1 7 5 8 3 4
5 4 9 6 5 1 4 4 6 4 4 4 4 4 4 4 8 0 3 5 9 3 7 7 4 9 9
2 6 2 5 2 0 4 4 4 4 4 4 4 4 4 4 5 3 7 0 4 1 1 3 1 3 4
0 9 9 3 1 3 4 4 4 4 4 4 4 6 4 4 6 2 8 8 0 8 9 2 0 5 5
2 6 3 9 2 0 4 4 4 4 4 4 4 4 4 4 5 5 1 6 9 7 3 1 9 3 5
6 2 0 2 4 6 4 4 4 4 4 6 4 4 4 4 7 4 1 1 1 7 0 9 0 5 1
9 1 9 8 5 8 4 4 4 4 4 4 4 4 4 4 2 6 3 2 5 7 1 2 8 8 6
4 9 6 7 2 7 4 4 4 4 4 4 4 4 4 4 7 8 9 8 5 5 7 7 6 3 7""".split('\n')

find_central_region(example1)