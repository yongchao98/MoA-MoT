def analyze_shift(input_str, output_str):
    input_list = [int(x) for x in input_str.split()]
    output_list = [int(x) for x in output_str.split()]
    
    # Find the non-zero sequences
    def find_sequence(arr):
        start = -1
        for i, val in enumerate(arr):
            if val != 0 and start == -1:
                start = i
            elif val == 0 and start != -1:
                return (start, i-1)
        return (start, len(arr)-1) if start != -1 else (-1, -1)
    
    in_start, in_end = find_sequence(input_list)
    out_start, out_end = find_sequence(output_list)
    
    if in_start != -1 and out_start != -1:
        shift = out_start - in_start
        length = out_end - out_start
        print(f"Shift: {shift}, Length: {length+1}")

# Test examples
print("Example 1:")
analyze_shift("0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0",
              "0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0")

print("\nExample 2:")
analyze_shift("0 0 0 0 0 0 0 0 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 0 0 0 0 0 0",
              "0 0 0 0 0 0 0 0 0 0 0 0 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 0 0")

print("\nExample 3:")
analyze_shift("0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0",
              "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0")