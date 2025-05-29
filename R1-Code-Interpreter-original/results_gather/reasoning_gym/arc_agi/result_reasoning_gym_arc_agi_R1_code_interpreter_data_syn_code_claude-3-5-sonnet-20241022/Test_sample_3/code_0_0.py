def analyze_pattern(input_examples, output_examples):
    # Convert examples to more readable format
    for i in range(len(input_examples)):
        print(f"\nExample {i+1}:")
        input_grid = input_examples[i]
        output_grid = output_examples[i]
        
        # Print positions of non-6 numbers in input
        print("Non-6 positions in input:")
        for row in range(6):
            for col in range(6):
                if input_grid[row][col] != 6:
                    print(f"Value {input_grid[row][col]} at position ({row}, {col})")
        
        # Print output grid
        print("\nOutput grid:")
        for row in output_grid:
            print(row)

# Define example grids
input_examples = [
    [[6,0,6,9,6,9],
     [6,6,6,6,6,6],
     [6,6,6,1,6,9],
     [6,6,6,6,6,6],
     [6,6,6,6,6,1],
     [6,6,6,6,6,6]],
    
    [[6,6,6,2,6,0],
     [6,6,6,6,6,6],
     [6,1,6,6,6,6],
     [6,6,6,6,6,6],
     [6,6,6,2,6,9],
     [6,6,6,6,6,6]],
    
    [[6,3,6,6,6,1],
     [6,6,6,6,6,6],
     [6,6,6,9,6,6],
     [6,6,6,6,6,6],
     [6,3,6,6,6,6],
     [6,6,6,6,6,6]]
]

output_examples = [
    [[0,9,9],
     [6,1,9],
     [6,6,1]],
    
    [[6,2,0],
     [1,6,6],
     [6,2,9]],
    
    [[3,6,1],
     [6,9,6],
     [3,6,6]]
]

analyze_pattern(input_examples, output_examples)