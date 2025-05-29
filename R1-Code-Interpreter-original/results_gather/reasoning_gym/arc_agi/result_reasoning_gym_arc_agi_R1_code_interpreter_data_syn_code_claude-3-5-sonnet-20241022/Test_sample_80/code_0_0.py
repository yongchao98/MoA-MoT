def analyze_pattern(input_grid, output_grid):
    # Convert input and output to 2D arrays
    input_arr = [list(map(int, row.split())) for row in input_grid.split('\n') if row]
    output_arr = [list(map(int, row.split())) for row in output_grid.split('\n') if row]
    
    # Find positions of 9's in output
    nine_positions = []
    for i in range(len(output_arr)):
        for j in range(len(output_arr[i])):
            if output_arr[i][j] == 9:
                nine_positions.append((i, j))
    
    # Count different numbers in input grid
    counts = {}
    for row in input_arr:
        for num in row:
            counts[num] = counts.get(num, 0) + 1
            
    return nine_positions, counts

# Example inputs and outputs
example1_in = """0 0 3 0 3 3 3
0 0 3 3 3 3 3
3 3 3 3 3 2 2
3 3 0 0 3 2 2
3 3 0 0 3 3 3"""

example1_out = """3 3 3
3 3 3
3 3 9"""

example2_in = """3 3 3 0 0
0 0 3 0 0
0 0 3 3 3
3 3 2 2 3
2 3 2 2 3"""

example2_out = """3 3 3
3 3 3
3 3 9"""

# Analyze patterns
pos1, counts1 = analyze_pattern(example1_in, example1_out)
pos2, counts2 = analyze_pattern(example2_in, example2_out)

print("Example 1 nine positions:", pos1)
print("Example 1 counts:", counts1)
print("\nExample 2 nine positions:", pos2)
print("Example 2 counts:", counts2)