def find_groups(input_seq, output_seq):
    n = len(input_seq)
    for group_size in range(1, n+1):
        # Try different group sizes
        if n % group_size == 0:
            groups_in = [input_seq[i:i+group_size] for i in range(0, n, group_size)]
            groups_out = [output_seq[i:i+group_size] for i in range(0, n, group_size)]
            print(f"\nGroup size {group_size}:")
            print(f"Input groups: {groups_in}")
            print(f"Output groups: {groups_out}")

# Test with Example 1
input_seq = [9, 6, 2, 6, 2, 5, 5, 8, 6, 8, 5, 7, 1, 7, 8]
output_seq = [8, 7, 1, 7, 5, 8, 6, 8, 5, 5, 2, 6, 2, 6, 9]
print("Analysis of Example 1:")
find_groups(input_seq, output_seq)