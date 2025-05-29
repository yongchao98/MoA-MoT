def analyze_pattern(input_str, output_str):
    input_nums = [int(x) for x in input_str.split()]
    output_nums = [int(x) for x in output_str.split()]
    
    # Find non-zero sequences
    first = next(i for i, x in enumerate(input_nums) if x != 0)
    last = len(input_nums) - next(i for i, x in enumerate(reversed(input_nums)) if x != 0) - 1
    
    input_seq = input_nums[first:last+1]
    output_seq = output_nums[first:last+1]
    
    # Split into chunks of 8
    def split_chunks(seq):
        return [seq[i:i+8] for i in range(0, len(seq), 8)]
    
    input_chunks = split_chunks(input_seq)
    output_chunks = split_chunks(output_seq)
    
    print("Input chunks:", input_chunks)
    print("Output chunks:", output_chunks)
    
    # Analyze each chunk transformation
    for i, (in_chunk, out_chunk) in enumerate(zip(input_chunks, output_chunks)):
        print(f"\nChunk {i+1}:")
        print(f"Input:  {in_chunk}")
        print(f"Output: {out_chunk}")

# Analyze examples
example1_in = "0 0 0 0 0 0 0 0 0 0 0 9 6 2 6 2 5 5 8 6 8 5 7 1 7 8 0 0"
example1_out = "0 0 0 0 0 0 0 0 0 0 0 8 7 1 7 5 8 6 8 5 5 2 6 2 6 9 0 0"
print("Example 1 analysis:")
analyze_pattern(example1_in, example1_out)