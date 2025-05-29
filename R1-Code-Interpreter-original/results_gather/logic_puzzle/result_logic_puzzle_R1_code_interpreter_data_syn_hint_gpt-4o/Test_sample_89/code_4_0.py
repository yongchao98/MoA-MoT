from itertools import combinations

grid = [
    [8, 12, 9, 3, -11, 9],
    [12, -12, -4, 11, 10, 7],
    [-8, 6, 9, 10, 1, 3],
    [5, 11, 6, 7, -9, 11],
    [6, 1, 11, 12, 4, -3],
    [2, 9, 4, 3, 10, 9]
]

# Start with the given number at position (0,5)
selected_positions = [(0, 5)]

# Function to check if the current selection is valid
def is_valid_selection(selected_positions):
    row_negatives = [0] * 6
    col_negatives = [0] * 6
    
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_negatives[r] += 1
            col_negatives[c] += 1
    
    for r in range(6):
        if row_negatives[r] % 2 != 0:
            return False
    for c in range(6):
        if col_negatives[c] % 2 != 0:
            return False
    return True

# Try all combinations of 11 more numbers
all_positions = [(r, c) for r in range(6) for c in range(6) if (r, c) != (0, 5)]
for combination in combinations(all_positions, 11):
    current_selection = selected_positions + list(combination)
    if is_valid_selection(current_selection):
        selected_positions = current_selection
        break

# Output the selected positions
print(', '.join(f"{r} {c}" for r, c in selected_positions))