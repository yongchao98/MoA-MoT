def analyze_transformation(input_str, output_str):
    input_grid = [int(x) for x in input_str.split()]
    output_grid = [int(x) for x in output_str.split()]
    
    # Find first non-zero sequence and its position
    def find_nonzero_sequence(grid):
        start = None
        length = 0
        for i, num in enumerate(grid):
            if num != 0:
                if start is None:
                    start = i
                length += 1
            elif start is not None and num == 0:
                return start, length
        return start, length if start is not None else 0

    in_start, in_length = find_nonzero_sequence(input_grid)
    out_start, out_length = find_nonzero_sequence(output_grid)
    
    print(f"Input sequence starts at {in_start}, length {in_length}")
    print(f"Output sequence starts at {out_start}, length {out_length}")
    print("Input sequence:", input_grid[in_start:in_start+in_length])
    print("Output sequence:", output_grid[out_start:out_start+out_length])
    print("Shift amount:", out_start - in_start if in_start is not None and out_start is not None else "N/A")
    print("-" * 50)

# Analyze examples
print("Example 1:")
analyze_transformation("0 0 0 0 0 0 0 0 0 0 0 0 0 0 9 3 1 4 0 0 0 0 0 0 0 0 0 0",
                      "0 0 0 0 0 0 0 0 0 0 0 9 3 1 4 0 0 0 0 0 0 0 0 0 0 0 0 0")

print("Example 2:")
analyze_transformation("2 6 6 9 2 5 9 0 0 0 0 0 0 0 0 0 0 0 8 7 9 6 5 9 2 7 8 2",
                      "9 2 5 9 0 0 0 0 0 0 0 0 0 0 0 8 7 9 6 5 9 2 7 8 2 2 6 6")

print("Example 3:")
analyze_transformation("0 0 0 0 0 6 5 3 2 2 5 6 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
                      "0 0 6 5 3 2 2 5 6 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0")

# Test input
test_input = "6 8 4 9 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9 2"
print("\nTest input analysis:")
analyze_transformation(test_input, "0 " * 28)  # dummy output for analysis