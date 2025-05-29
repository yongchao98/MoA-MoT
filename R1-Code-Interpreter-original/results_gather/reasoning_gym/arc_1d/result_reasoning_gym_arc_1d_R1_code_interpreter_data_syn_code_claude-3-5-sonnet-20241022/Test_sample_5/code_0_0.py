def analyze_pattern(input_grid, output_grid):
    # Convert string to list of integers
    input_nums = [int(x) for x in input_grid.split()]
    output_nums = [int(x) for x in output_grid.split()]
    
    # Find the non-zero sequences
    def get_nonzero_seq(arr):
        non_zero = [x for x in arr if x != 0]
        start = next(i for i, x in enumerate(arr) if x != 0)
        return non_zero, start
    
    in_nonzero, in_start = get_nonzero_seq(input_nums)
    out_nonzero, out_start = get_nonzero_seq(output_nums)
    
    print(f"Input sequence: starts at {in_start}, length {len(in_nonzero)}")
    print(f"Output sequence: starts at {out_start}, length {len(out_nonzero)}")

# Test examples
examples = [
    ("0 0 0 0 0 0 0 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6",
     "0 0 0 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 0 0 0 0"),
    ("0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0",
     "0 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 0"),
    ("0 0 0 0 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 0 0 0",
     "5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 0 0 0 0 0 0 0")
]

for i, (inp, out) in enumerate(examples, 1):
    print(f"\nExample {i}:")
    analyze_pattern(inp, out)