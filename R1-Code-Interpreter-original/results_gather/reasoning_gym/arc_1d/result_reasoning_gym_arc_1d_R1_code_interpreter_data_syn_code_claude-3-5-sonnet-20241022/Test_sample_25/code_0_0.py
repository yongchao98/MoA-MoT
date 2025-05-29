def analyze_shift(input_str, output_str):
    input_list = [int(x) for x in input_str.split()]
    output_list = [int(x) for x in output_str.split()]
    
    # Find non-zero sequences
    def get_nonzero_seq(lst):
        start = next((i for i, x in enumerate(lst) if x != 0), None)
        if start is None: return None, None, []
        end = len(lst) - next((i for i, x in enumerate(reversed(lst)) if x != 0), 0)
        return start, end, lst[start:end]
    
    in_start, in_end, in_seq = get_nonzero_seq(input_list)
    out_start, out_end, out_seq = get_nonzero_seq(output_list)
    
    # Compare sequences and determine shift
    if in_seq == out_seq:
        shift = out_start - in_start
        print(f"Sequence: {in_seq}")
        print(f"Input position: {in_start}")
        print(f"Output position: {out_start}")
        print(f"Shift: {shift}")

# Analyze examples
print("Example 1:")
analyze_shift("0 0 0 0 0 0 0 9 0 0 0 0 0 0 0 0 0 0 0 0 0",
             "0 0 0 0 0 0 0 0 0 9 0 0 0 0 0 0 0 0 0 0 0")
print("\nExample 2:")
analyze_shift("0 0 0 0 0 0 0 5 7 7 1 3 9 6 7 7 2 8 0 0 0",
             "0 0 0 0 0 0 0 0 0 5 7 7 1 3 9 6 7 7 2 8 0")
print("\nExample 3:")
analyze_shift("0 3 3 3 7 9 3 7 9 3 7 6 9 9 4 4 5 7 4 0 0",
             "0 0 0 3 3 3 7 9 3 7 9 3 7 6 9 9 4 4 5 7 4")