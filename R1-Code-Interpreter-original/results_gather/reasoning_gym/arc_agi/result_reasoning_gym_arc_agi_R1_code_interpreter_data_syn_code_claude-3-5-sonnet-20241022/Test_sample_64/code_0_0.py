def analyze_transformation(input_grid, output_grid):
    print("Analyzing transformation:")
    print("Input grid:")
    for row in input_grid:
        print(row)
    print("\nOutput grid:")
    for row in output_grid:
        print(row)
    
    # Count occurrences of each number in input and output
    input_counts = {}
    output_counts = {}
    for i in range(3):
        for j in range(3):
            input_counts[input_grid[i][j]] = input_counts.get(input_grid[i][j], 0) + 1
            output_counts[output_grid[i][j]] = output_counts.get(output_grid[i][j], 0) + 1
    
    print("\nNumber counts (should be same in input and output):")
    print("Input counts:", input_counts)
    print("Output counts:", output_counts)

# Example 1
input1 = [[4,9,9], [4,4,9], [9,4,9]]
output1 = [[9,4,4], [4,4,9], [9,9,9]]

# Example 2
input2 = [[3,3,3], [4,4,3], [4,4,3]]
output2 = [[4,4,3], [4,4,3], [3,3,3]]

analyze_transformation(input1, output1)
print("\n-------------------\n")
analyze_transformation(input2, output2)