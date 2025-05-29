def analyze_pattern(input_grid, output_grid):
    # Find the first and last positions of non-zero numbers in both grids
    in_start = next(i for i, x in enumerate(input_grid) if x != 0)
    in_end = len(input_grid) - next(i for i, x in enumerate(reversed(input_grid)) if x != 0) - 1
    out_start = next(i for i, x in enumerate(output_grid) if x != 0)
    
    # Print analysis
    block_length = in_end - in_start + 1
    shift = in_start - out_start
    print(f"Block length: {block_length}")
    print(f"Original position: {in_start}-{in_end}")
    print(f"New position: {out_start}-{out_start+block_length-1}")
    print(f"Shift left by: {shift}")
    return block_length, shift

# Example inputs and outputs
ex1_in = [0,0,0,0,0,0,0,0,0,2,2,2,0,0,0,0]
ex1_out = [0,0,0,0,0,2,2,2,0,0,0,0,0,0,0,0]
ex3_in = [0,0,0,0,0,8,8,8,8,8,8,0,0,0,0,0]
ex3_out = [0,8,8,8,8,8,8,0,0,0,0,0,0,0,0,0]

print("Example 1 analysis:")
analyze_pattern(ex1_in, ex1_out)
print("\nExample 3 analysis:")
analyze_pattern(ex3_in, ex3_out)

# Test input
test_in = [0,0,0,0,7,7,7,7,7,7,7,7,7,0,0,0]

# Find the block in test input
test_start = next(i for i, x in enumerate(test_in) if x != 0)
test_end = len(test_in) - next(i for i, x in enumerate(reversed(test_in)) if x != 0) - 1
block_length = test_end - test_start + 1

# Generate output based on pattern
test_out = [0] * 16
new_start = 1  # Based on observed pattern
for i in range(block_length):
    test_out[new_start + i] = 7

print("\nPredicted output:")
print(" ".join(map(str, test_out)))