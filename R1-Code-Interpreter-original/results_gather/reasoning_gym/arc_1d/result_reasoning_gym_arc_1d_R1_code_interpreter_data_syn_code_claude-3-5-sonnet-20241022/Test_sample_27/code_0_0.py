def analyze_pattern(input_seq, output_seq):
    print("Length comparison:", len(input_seq) == len(output_seq))
    
    # Check if output is a shifted version of input
    for i in range(len(input_seq)):
        shifted = input_seq[i:] + input_seq[:i]
        if shifted == output_seq:
            print(f"Found shift pattern: shift left by {i} positions")
            return i

# Test examples
example1_in = [1,1,5,3,1,0,3,5,9,8,7,8]
example1_out = [1,5,3,1,0,3,5,9,8,7,8,1]

example2_in = [9,3,6,4,1,2,5,6,1,6,2,2]
example2_out = [3,6,4,1,2,5,6,1,6,2,2,9]

example3_in = [8,2,1,7,3,6,2,2,5,5,2,5]
example3_out = [2,1,7,3,6,2,2,5,5,2,5,8]

print("Example 1:")
analyze_pattern(example1_in, example1_out)
print("\nExample 2:")
analyze_pattern(example2_in, example2_out)
print("\nExample 3:")
analyze_pattern(example3_in, example3_out)