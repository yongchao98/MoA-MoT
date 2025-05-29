def analyze_transformation(input_str, output_str):
    # Convert strings to lists of numbers
    input_nums = [int(x) for x in input_str.split()]
    output_nums = [int(x) for x in output_str.split()]
    
    # Find the first and last non-zero positions
    first_nonzero = next(i for i, x in enumerate(input_nums) if x != 0)
    last_nonzero = len(input_nums) - next(i for i, x in enumerate(reversed(input_nums)) if x != 0) - 1
    
    # Extract the non-zero sequences
    input_seq = input_nums[first_nonzero:last_nonzero+1]
    output_seq = output_nums[first_nonzero:last_nonzero+1]
    
    print(f"Input sequence: {input_seq}")
    print(f"Output sequence: {output_seq}")
    print(f"Length: {len(input_seq)}")
    print("---")

# Analyze each example
example1_in = "0 0 0 0 0 0 0 0 0 0 0 9 6 2 6 2 5 5 8 6 8 5 7 1 7 8 0 0"
example1_out = "0 0 0 0 0 0 0 0 0 0 0 8 7 1 7 5 8 6 8 5 5 2 6 2 6 9 0 0"
example2_in = "0 3 5 6 2 8 5 1 4 6 7 1 1 1 4 4 2 1 1 8 5 1 8 4 0 0 0 0"
example2_out = "0 4 8 1 5 8 1 1 2 4 4 1 1 1 7 6 4 1 5 8 2 6 5 3 0 0 0 0"
example3_in = "0 0 0 3 3 5 6 2 7 2 4 6 6 6 1 8 0 0 0 0 0 0 0 0 0 0 0 0"
example3_out = "0 0 0 8 1 6 6 6 4 2 7 2 6 5 3 3 0 0 0 0 0 0 0 0 0 0 0 0"

print("Example 1:")
analyze_transformation(example1_in, example1_out)
print("Example 2:")
analyze_transformation(example2_in, example2_out)
print("Example 3:")
analyze_transformation(example3_in, example3_out)