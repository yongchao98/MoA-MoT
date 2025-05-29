def analyze_test_input(input_grid):
    # Find positions of non-zero numbers
    non_zeros = [(i, x) for i, x in enumerate(input_grid) if x != 0]
    print("Non-zero positions:", non_zeros)
    
    # Count total numbers
    total_len = len(input_grid)
    zeros_count = input_grid.count(0)
    print(f"Total length: {total_len}")
    print(f"Zeros count: {zeros_count}")
    
    # Analyze gaps
    if non_zeros:
        before_gap = non_zeros[0][0]
        after_gap = total_len - non_zeros[-1][0] - 1
        print(f"Gap before first non-zero: {before_gap}")
        print(f"Gap after last non-zero: {after_gap}")

test_input = [0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0]
analyze_test_input(test_input)