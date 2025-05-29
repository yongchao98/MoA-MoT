def analyze_pattern(input_grid, output_grid):
    base_num = input_grid[0]
    print(f"\nAnalyzing pattern for number {base_num}:")
    
    # Find positions of non-zero numbers
    input_positions = [i for i, x in enumerate(input_grid) if x != 0]
    output_positions = [i for i, x in enumerate(output_grid) if x != 0]
    
    print("Input non-zero positions:", input_positions)
    print("Output non-zero positions:", output_positions)
    
    # Print as 5x5 grid
    def print_as_grid(positions, size=5):
        grid = [[0 for _ in range(size)] for _ in range(size)]
        for pos in positions:
            if pos < size * size:
                grid[pos // size][pos % size] = 1
        for row in grid:
            print(row)
    
    print("\nInput grid:")
    print_as_grid(input_positions)
    print("\nOutput grid:")
    print_as_grid(output_positions)

# Example 1
input1 = [2,2,2,2,2,0,0,0,2,0,0,0,0,0,0,2,0,0,0,0,0,0,0]
output1 = [2,2,2,2,2,0,2,2,2,2,2,0,0,2,2,2,2,2,0,0,0,0,0]

# Example 2
input2 = [3,3,3,3,3,0,0,0,0,3,0,0,0,0,0,0,0,3,0,0,0,0,0]
output2 = [3,3,3,3,3,0,0,3,3,3,3,3,0,0,0,3,3,3,3,3,0,0,0]

analyze_pattern(input1, output1)
analyze_pattern(input2, output2)