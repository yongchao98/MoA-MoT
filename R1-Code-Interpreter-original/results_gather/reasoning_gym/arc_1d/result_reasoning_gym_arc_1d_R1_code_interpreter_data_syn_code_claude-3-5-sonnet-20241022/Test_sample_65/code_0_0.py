def analyze_pattern(input_str, output_str):
    # Convert strings to lists of integers
    input_list = [int(x) for x in input_str.split()]
    output_list = [int(x) for x in output_str.split()]
    
    # Print analysis
    print(f"Input:  {input_list}")
    print(f"Output: {output_list}")
    print(f"Difference positions:", [i for i in range(len(input_list)) if input_list[i] != output_list[i]])
    print("---")

# Example cases
example1_in = "0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0"
example1_out = "0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0"

example2_in = "0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0"
example2_out = "0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0"

example3_in = "0 0 0 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0"
example3_out = "0 0 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0"

# Analyze examples
analyze_pattern(example1_in, example1_out)
analyze_pattern(example2_in, example2_out)
analyze_pattern(example3_in, example3_out)

# Function to apply the rule
def apply_rule(input_str):
    input_list = [int(x) for x in input_str.split()]
    output_list = input_list.copy()
    
    # Find the first non-zero number
    first_nonzero = -1
    for i in range(len(input_list)):
        if input_list[i] != 0:
            first_nonzero = i
            break
    
    if first_nonzero != -1:
        # Shift the sequence one position left
        for i in range(first_nonzero, len(input_list)-1):
            output_list[i] = input_list[i+1]
        output_list[-1] = 0  # Last position becomes 0
        
    return ' '.join(map(str, output_list))

# Test input
test_input = "8 8 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
print("Test output:", apply_rule(test_input))