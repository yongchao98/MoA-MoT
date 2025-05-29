def analyze_chunks_detailed(input_str, output_str):
    input_nums = [int(x) for x in input_str.split()]
    output_nums = [int(x) for x in output_str.split()]
    
    # Find non-zero sequences
    first = next(i for i, x in enumerate(input_nums) if x != 0)
    last = len(input_nums) - next(i for i, x in enumerate(reversed(input_nums)) if x != 0) - 1
    
    input_seq = input_nums[first:last+1]
    output_seq = output_nums[first:last+1]
    
    # Split into chunks of 8
    input_chunks = [input_seq[i:i+8] for i in range(0, len(input_seq), 8)]
    output_chunks = [output_seq[i:i+8] for i in range(0, len(output_seq), 8)]
    
    for i, (in_chunk, out_chunk) in enumerate(zip(input_chunks, output_chunks)):
        print(f"\nChunk {i+1}:")
        print(f"Input:  {in_chunk}")
        print(f"Output: {out_chunk}")
        # Check if it's a reverse
        print(f"Is reverse: {in_chunk[::-1] == out_chunk}")
        # Check if numbers are preserved
        print(f"Input sorted:  {sorted(in_chunk)}")
        print(f"Output sorted: {sorted(out_chunk)}")

# Analyze Example 1
example1_in = "0 0 0 0 0 0 0 0 0 0 0 9 6 2 6 2 5 5 8 6 8 5 7 1 7 8 0 0"
example1_out = "0 0 0 0 0 0 0 0 0 0 0 8 7 1 7 5 8 6 8 5 5 2 6 2 6 9 0 0"
print("Example 1 detailed analysis:")
analyze_chunks_detailed(example1_in, example1_out)