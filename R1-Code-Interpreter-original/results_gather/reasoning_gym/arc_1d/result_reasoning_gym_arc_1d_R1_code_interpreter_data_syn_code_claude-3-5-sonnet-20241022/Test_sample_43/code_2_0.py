def analyze_positions(input_str, output_str):
    input_nums = [int(x) for x in input_str.split()]
    output_nums = [int(x) for x in output_str.split()]
    
    # Find non-zero sequences
    first = next(i for i, x in enumerate(input_nums) if x != 0)
    last = len(input_nums) - next(i for i, x in enumerate(reversed(input_nums)) if x != 0) - 1
    
    input_seq = input_nums[first:last+1]
    output_seq = output_nums[first:last+1]
    
    # Create position mapping
    positions = []
    for num in input_seq:
        # Find where each input number appears in output
        out_pos = output_seq.index(num)
        positions.append(out_pos)
    
    print(f"Input sequence: {input_seq}")
    print(f"Output sequence: {output_seq}")
    print(f"Position mapping: {positions}")

# Analyze Example 1
example1_in = "0 0 0 0 0 0 0 0 0 0 0 9 6 2 6 2 5 5 8 6 8 5 7 1 7 8 0 0"
example1_out = "0 0 0 0 0 0 0 0 0 0 0 8 7 1 7 5 8 6 8 5 5 2 6 2 6 9 0 0"
print("Example 1 analysis:")
analyze_positions(example1_in, example1_out)