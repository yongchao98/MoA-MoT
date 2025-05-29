def analyze_pattern(input_grid, output_grid):
    # Find non-zero sequences
    input_nonzero = [(i, val) for i, val in enumerate(input_grid) if val != 0]
    output_nonzero = [(i, val) for i, val in enumerate(output_grid) if val != 0]
    
    # Print analysis
    print(f"Input sequence: {input_nonzero}")
    print(f"Output sequence: {output_nonzero}")
    print(f"Shift: {output_nonzero[0][0] - input_nonzero[0][0] if input_nonzero and output_nonzero else 'N/A'}")
    print("---")

# Example cases
ex1_in = [7,7,7,7,7,7,7,7,7,7,7,7,0,0,0,0]
ex1_out = [0,0,0,7,7,7,7,7,7,7,7,7,7,7,7,0]

ex2_in = [3,3,3,3,3,3,0,0,0,0,0,0,0,0,0,0]
ex2_out = [0,0,0,3,3,3,3,3,3,0,0,0,0,0,0,0]

ex3_in = [0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0]
ex3_out = [0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0]

test_in = [0,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0]

print("Pattern Analysis:")
analyze_pattern(ex1_in, ex1_out)
analyze_pattern(ex2_in, ex2_out)
analyze_pattern(ex3_in, ex3_out)
analyze_pattern(test_in, [0]*16)  # For visualization