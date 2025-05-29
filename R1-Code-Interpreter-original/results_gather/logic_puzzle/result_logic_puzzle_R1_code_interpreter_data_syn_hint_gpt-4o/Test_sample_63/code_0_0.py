from itertools import combinations

# The given grid
grid = [
    [8, -3, 5, 1, 3],
    [8, 5, -6, 3, 4],
    [5, -8, 7, -1, -8],
    [-2, -3, -4, -7, 1],
    [-1, 8, -8, -4, 8]
]

# The initial selected position
initial_position = (0, 4)

# Function to check if the product of numbers in each row and column is positive
def is_valid_selection(selected_positions):
    # Create a 5x5 matrix to mark selected positions
    selected_matrix = [[0] * 5 for _ in range(5)]
    for r, c in selected_positions:
        selected_matrix[r][c] = 1
    
    # Check rows
    for r in range(5):
        product = 1
        for c in range(5):
            if selected_matrix[r][c] == 1:
                product *= grid[r][c]
        if product <= 0:
            return False
    
    # Check columns
    for c in range(5):
        product = 1
        for r in range(5):
            if selected_matrix[r][c] == 1:
                product *= grid[r][c]
        if product <= 0:
            return False
    
    return True

# Generate all possible combinations of 8 positions (since one is already selected)
all_positions = [(r, c) for r in range(5) for c in range(5) if (r, c) != initial_position]
for combination in combinations(all_positions, 8):
    selected_positions = [initial_position] + list(combination)
    if is_valid_selection(selected_positions):
        # Format the output as required
        result = ', '.join(f"{r} {c}" for r, c in selected_positions)
        print(f"<<<{result}>>>")
        break