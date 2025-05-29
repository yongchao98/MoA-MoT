import itertools

# The grid
grid = [
    [-1, -6, 5, 5, -8],
    [1, -5, 5, -1, -4],
    [5, 8, 6, 5, 3],
    [7, 5, 2, 3, 5],
    [6, 4, 6, 7, 3]
]

# Pre-select the given number
selected_positions = [(0, 4)]

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(selected_positions):
    # Create a 5x5 grid to mark selected positions
    selected_grid = [[0] * 5 for _ in range(5)]
    for r, c in selected_positions:
        selected_grid[r][c] = 1
    
    # Check rows
    for r in range(5):
        product = 1
        for c in range(5):
            if selected_grid[r][c] == 1:
                product *= grid[r][c]
        if product <= 0:
            return False
    
    # Check columns
    for c in range(5):
        product = 1
        for r in range(5):
            if selected_grid[r][c] == 1:
                product *= grid[r][c]
        if product <= 0:
            return False
    
    return True

# Try to find a valid selection
def find_valid_selection():
    # All possible positions except the pre-selected one
    all_positions = [(r, c) for r in range(5) for c in range(5) if (r, c) != (0, 4)]
    
    # Try all combinations of 8 additional positions
    for additional_positions in itertools.combinations(all_positions, 8):
        current_selection = selected_positions + list(additional_positions)
        if is_valid_selection(current_selection):
            return current_selection

# Find and print the valid selection
valid_selection = find_valid_selection()
print(valid_selection)