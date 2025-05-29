import itertools

grid = [
    [5, 1, 10, 11, 1, 1],
    [7, 9, -1, 12, -9, 8],
    [3, 10, 1, -1, -11, 1],
    [3, -10, -4, -9, 3, 6],
    [3, -2, 3, -10, 2, 3],
    [-9, -3, 8, 4, 6, 6]
]

def is_valid_selection(selected_positions):
    row_neg_count = [0] * 6
    col_neg_count = [0] * 6
    
    for r, c in selected_positions:
        if grid[r][c] < 0:
            row_neg_count[r] += 1
            col_neg_count[c] += 1
    
    for count in row_neg_count + col_neg_count:
        if count % 2 != 0:
            return False
    return True

def find_selection():
    initial_position = (0, 5)
    all_positions = [(i, j) for i in range(6) for j in range(6) if (i, j) != initial_position]
    
    for combination in itertools.combinations(all_positions, 11):
        selected_positions = [initial_position] + list(combination)
        if is_valid_selection(selected_positions):
            return selected_positions

selected_positions = find_selection()
print(selected_positions)