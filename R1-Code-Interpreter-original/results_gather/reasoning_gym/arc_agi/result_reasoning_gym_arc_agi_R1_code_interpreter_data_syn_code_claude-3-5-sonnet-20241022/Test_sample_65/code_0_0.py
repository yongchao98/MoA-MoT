def analyze_pattern(input_grid, output_grid):
    # Find the special number in the last row
    last_row = input_grid[-1]
    special_pos = None
    special_num = None
    for i, num in enumerate(last_row):
        if num != 1:
            special_pos = i
            special_num = num
            break
    
    # Analyze the pattern in output
    n = len(output_grid)
    m = len(output_grid[0])
    
    print(f"Special number {special_num} at position {special_pos}")
    print("\nVertical patterns from special position:")
    for i in range(special_pos-2, special_pos+3):
        if 0 <= i < m:
            column = [output_grid[j][i] for j in range(n)]
            print(f"Column {i}: {column}")

# Example 1
input1 = [[1]*10 for _ in range(9)] + [[1,1,1,1,1,1,1,1,6,1]]
output1 = [
    [6,1,6,5,6,1,6,5,6,1],
    [6,1,6,1,6,1,6,1,6,1],
    [6,1,6,1,6,1,6,1,6,1],
    [6,1,6,1,6,1,6,1,6,1],
    [6,1,6,1,6,1,6,1,6,1],
    [6,1,6,1,6,1,6,1,6,1],
    [6,1,6,1,6,1,6,1,6,1],
    [6,1,6,1,6,1,6,1,6,1],
    [6,1,6,1,6,1,6,1,6,1],
    [6,5,6,1,6,5,6,1,6,1]
]

print("Analysis of Example 1:")
analyze_pattern(input1, output1)