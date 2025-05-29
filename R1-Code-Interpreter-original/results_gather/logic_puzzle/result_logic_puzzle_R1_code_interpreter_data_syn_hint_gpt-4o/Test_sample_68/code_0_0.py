grid = [
    [9, 7, 4, 5, 5],
    [-7, 4, -3, -4, 5],
    [7, 1, 3, 2, 8],
    [4, 9, 2, 1, 6],
    [-4, 6, -9, 7, 3]
]

# Start with the given number
selected_positions = [(0, 4)]
selected_numbers = [grid[0][4]]

# Function to check if the current selection is valid
def is_valid_selection(selected_positions):
    row_count = [0] * 5
    col_count = [0] * 5
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_count[r] += 1
            col_count[c] += 1
    # Check if each row and column has an even number of negative numbers
    for count in row_count + col_count:
        if count % 2 != 0:
            return False
    return True

# Try to select 10 numbers
for r in range(5):
    for c in range(5):
        if (r, c) not in selected_positions:
            selected_positions.append((r, c))
            selected_numbers.append(grid[r][c])
            if len(selected_positions) == 10:
                if is_valid_selection(selected_positions):
                    break
                else:
                    selected_positions.pop()
                    selected_numbers.pop()
    if len(selected_positions) == 10 and is_valid_selection(selected_positions):
        break

# Output the selected positions
print(selected_positions)