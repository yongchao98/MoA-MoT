def analyze_pattern(input_grid, output_grid):
    # Count occurrences of each number
    unique_nums = sorted(set(input_grid))
    counts_in = {n: input_grid.count(n) for n in unique_nums}
    counts_out = {n: output_grid.count(n) for n in unique_nums}
    
    # Analyze the pattern
    print(f"Input counts: {counts_in}")
    print(f"Output counts: {counts_out}")
    
    # Print the positions of non-zero numbers
    non_zero = [n for n in unique_nums if n != 0]
    for n in non_zero:
        in_pos = [i for i, x in enumerate(input_grid) if x == n]
        out_pos = [i for i, x in enumerate(output_grid) if x == n]
        print(f"\nNumber {n}:")
        print(f"Input positions: {in_pos}")
        print(f"Output positions: {out_pos}")

# Example 1
input1 = [0,0,0,0,0,4,4,4,4,4,4,4,4,4,4,4,4,4]
output1 = [0,0,4,4,4,4,4,4,4,4,4,4,4,4,4,0,0,0]

print("Analysis of Example 1:")
analyze_pattern(input1, output1)