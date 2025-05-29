def analyze_shift(input_grid, output_grid):
    # Find first non-zero position in input and output
    def first_nonzero(arr):
        for i, x in enumerate(arr):
            if x != 0:
                return i
        return -1
    
    input_start = first_nonzero(input_grid)
    output_start = first_nonzero(output_grid)
    
    # Get non-zero sequences
    input_seq = [x for x in input_grid[input_start:] if x != 0]
    output_seq = [x for x in output_grid[output_start:] if x != 0]
    
    print(f"Input starts at: {input_start}")
    print(f"Output starts at: {output_start}")
    print(f"Shift amount: {output_start - input_start}")
    print(f"Input sequence: {input_seq}")
    print(f"Output sequence: {output_seq}")

# Test examples
example1_in = [0,6,6,9,0,0,0,0,0,0,0,0,0]
example1_out = [0,0,6,6,9,0,0,0,0,0,0,0,0]
example2_in = [0,0,0,7,8,9,6,9,9,8,0,0,0]
example2_out = [0,0,0,0,7,8,9,6,9,9,8,0,0]
example3_in = [0,0,0,0,0,0,7,2,7,0,0,0,0]
example3_out = [0,0,0,0,0,0,0,7,2,7,0,0,0]

print("Example 1:")
analyze_shift(example1_in, example1_out)
print("\nExample 2:")
analyze_shift(example2_in, example2_out)
print("\nExample 3:")
analyze_shift(example3_in, example3_out)