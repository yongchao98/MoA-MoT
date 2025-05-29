def analyze_pattern(input_grid, output_grid):
    # Count non-zero numbers
    input_nonzero = [x for x in input_grid if x != 0]
    output_nonzero = [x for x in output_grid if x != 0]
    
    # Create expected output based on our hypothesis
    expected = input_nonzero + [0] * (len(input_grid) - len(input_nonzero))
    
    # Compare with actual output
    print(f"Input non-zero count: {len(input_nonzero)}")
    print(f"Output non-zero count: {len(output_nonzero)}")
    print(f"Expected: {expected}")
    print(f"Actual:   {output_grid}")
    print(f"Match: {expected == output_grid}")

# Test examples
example1_in = [0,6,0,0,6,6,6,6,6,6,0,0,6,6,6,0,0,6,6,0,6,6,6,0]
example1_out = [6,0,6,6,0,6,6,6,0,6,6,6,0,6,6,6,6,6,6,0,0,0,0,0]

example2_in = [6,6,0,0,0,0,6,6,6,6,0,0,0,0,0,6,6,6,0,0,0,0,0,0]
example2_out = [6,6,0,6,6,6,0,6,6,6,6,0,0,0,0,0,0,0,0,0,0,0,0,0]

print("Example 1:")
analyze_pattern(example1_in, example1_out)
print("\nExample 2:")
analyze_pattern(example2_in, example2_out)