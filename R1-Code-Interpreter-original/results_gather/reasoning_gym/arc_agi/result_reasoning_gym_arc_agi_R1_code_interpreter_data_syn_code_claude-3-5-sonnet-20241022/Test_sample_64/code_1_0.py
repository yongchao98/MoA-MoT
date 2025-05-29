def analyze_positions(input_grid, output_grid):
    print("Analyzing grid positions:")
    print("Input grid:")
    for row in input_grid:
        print(row)
    print("\nOutput grid:")
    for row in output_grid:
        print(row)
    
    # Analyze row and column patterns
    print("\nRow sums:")
    print("Input:", [sum(row) for row in input_grid])
    print("Output:", [sum(row) for row in output_grid])
    
    # Column sums
    input_cols = [sum(input_grid[i][j] for i in range(3)) for j in range(3)]
    output_cols = [sum(output_grid[i][j] for i in range(3)) for j in range(3)]
    print("\nColumn sums:")
    print("Input:", input_cols)
    print("Output:", output_cols)

# Example 3
input3 = [[6,6,4], [4,4,4], [6,6,6]]
output3 = [[6,4,6], [6,4,6], [6,4,4]]

# Example 4
input4 = [[3,3,3], [3,3,3], [3,4,4]]
output4 = [[3,3,3], [4,3,3], [4,3,3]]

analyze_positions(input3, output3)
print("\n-------------------\n")
analyze_positions(input4, output4)