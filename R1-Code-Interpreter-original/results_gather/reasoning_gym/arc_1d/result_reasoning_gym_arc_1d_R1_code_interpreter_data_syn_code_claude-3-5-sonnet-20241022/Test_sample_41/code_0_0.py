def analyze_transformation(input_str, output_str):
    input_list = input_str.split()
    output_list = output_str.split()
    
    print(f"Input:  {input_list}")
    print(f"Output: {output_list}")
    print("Changes:")
    
    # Find first different position
    for i in range(len(input_list)):
        if input_list[i] != output_list[i]:
            print(f"First difference at position {i}")
            print(f"Input subsequence: {input_list[i:i+3]}")
            print(f"Output subsequence: {output_list[i:i+3]}")
            break
    print("---")

# Analyze all examples
example1_in = "0 0 2 9 0 1 1 9 9 0 1 0 0"
example1_out = "0 2 9 0 1 1 9 9 0 1 0 0 0"
example2_in = "2 6 0 0 4 0 0 0 0 0 0 0 0"
example2_out = "2 6 0 4 0 0 0 0 0 0 0 0 0"
example3_in = "2 0 4 3 8 0 0 5 7 7 1 0 0"
example3_out = "2 4 3 8 0 0 5 7 7 1 0 0 0"

analyze_transformation(example1_in, example1_out)
analyze_transformation(example2_in, example2_out)
analyze_transformation(example3_in, example3_out)