def analyze_grid(input_grid):
    # Find the column with 1s (divider)
    divider_col = None
    for j in range(len(input_grid[0])):
        if input_grid[0][j] == '1':
            divider_col = j
            break
    
    # Find rows with special numbers (4 or 9) and their positions
    special_sequences = []
    for i in range(len(input_grid)):
        row = input_grid[i]
        for num in ['4', '9']:
            # Check left side of divider
            left_sequence = ''.join(row[:divider_col]).count(num * 3)
            # Check right side of divider
            right_sequence = ''.join(row[divider_col+1:]).count(num * 2)
            if left_sequence > 0 or right_sequence > 0:
                special_sequences.append((i, num))
    
    if not special_sequences:
        return None
    
    # Sort by row number to find the bottommost special sequence
    special_sequences.sort(key=lambda x: x[0], reverse=True)
    # The bottommost special sequence determines the output
    print(f"Special sequences found (row, number): {special_sequences}")
    return special_sequences[0][1]

# Test input
test_input = [
    "8 8 8 8 8 1 8 8 8 8 8".split(),
    "8 8 8 8 8 1 8 8 8 8 8".split(),
    "8 9 9 9 9 1 8 8 8 9 9".split(),
    "8 8 8 8 8 1 8 8 8 8 8".split(),
    "8 8 8 8 8 1 8 8 8 8 8".split(),
    "8 8 8 8 8 1 8 8 8 8 8".split(),
    "8 4 4 4 4 1 8 8 4 4 4".split(),
    "8 8 8 8 8 1 8 8 8 8 8".split(),
    "8 8 8 8 8 1 8 8 8 8 8".split()
]

result = analyze_grid(test_input)