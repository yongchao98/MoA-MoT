def analyze_all_examples():
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
    
    # Example 2
    input2 = [[1]*10 for _ in range(9)] + [[1,1,1,1,4,1,1,1,1,1]]
    output2 = [
        [4,1,4,5,4,1,1,1,1,1],
        [4,1,4,1,4,1,1,1,1,1],
        [4,1,4,1,4,1,1,1,1,1],
        [4,1,4,1,4,1,1,1,1,1],
        [4,1,4,1,4,1,1,1,1,1],
        [4,1,4,1,4,1,1,1,1,1],
        [4,1,4,1,4,1,1,1,1,1],
        [4,1,4,1,4,1,1,1,1,1],
        [4,1,4,1,4,1,1,1,1,1],
        [4,5,4,1,4,1,1,1,1,1]
    ]
    
    for example_num, (input_grid, output_grid) in enumerate([
        (input1, output1), (input2, output2)
    ], 1):
        print(f"\nExample {example_num}:")
        # Find special number and its position
        special_pos = next(i for i, x in enumerate(input_grid[-1]) if x != 1)
        special_num = input_grid[-1][special_pos]
        print(f"Special number {special_num} at position {special_pos}")
        
        # Print the first row pattern
        print("First row:", output_grid[0])
        # Print the last row pattern
        print("Last row:", output_grid[-1])
        # Print positions of 5s
        five_positions = [(i, j) for i in range(10) for j in range(10) if output_grid[i][j] == 5]
        print("Positions of 5:", five_positions)

analyze_all_examples()